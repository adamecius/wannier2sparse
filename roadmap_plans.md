# wannier2sparse — Planes de ejecución completos (autocontenido)

Documento de referencia único para el desarrollo full-C++ de
`adamecius/wannier2sparse`. Pensado para sobrevivir a una compactación de
contexto: incluye estado del repo, arquitectura, invariantes y los 8 planes
acotados en detalle.

Repo verificado contra `master` = `c6656ec` (13-jun-2026). Último cambio:
solo se añadió `.github/workflows/ci.yml`; ningún fuente tocado.

---

## 0. Contexto y objetivo

wannier2sparse es una herramienta C++11 (Eigen header-only, CMake ≥ 3.14, parte
histórica de LinQT) que expande un modelo tight-binding de Wannier90 a un
supercell y exporta el Hamiltoniano y operadores como matrices sparse CSR, listas
para KPM/Chebyshev.

**Objetivo**: convertirla en el generador de operadores en espacio real que
produce los **goldens** del refactor de LinQT, manteniéndola **full C++** (un solo
toolchain, ruta byte-estable en un sitio). Los wrappers Python se dejan para
después. wannier2sparse **desaparecerá de dentro de LinQT**; esta versión
standalone es la que se usa.

**Operadores objetivo finales**: velocidad/corriente, corriente de spin, operador
de spin, momento angular orbital.

### Estado actual del código (mapa)

```text
src/wannier_parser.cpp     read_wannier_file, read_xyz_file, read_unit_cell_file
src/hopping_list.cpp       create_hopping_list, wrap_in_supercell, save_supercell_as_csr
src/main.cpp               CLI: build_operator() despacha a tbmodel::*; expande y exporta
src/sparse_matrix.cpp
include/hopping_list.hpp    struct hopping_list, save_hopping_list_as_csr (inline)
include/tbmodel.hpp         class tbmodel: builders de operadores
include/wannier_parser.hpp
include/sparse_matrix.hpp   SparseMatrix_t = Eigen::SparseMatrix<complex<double>,RowMajor>
include/w2sp_arguments.hpp  class W2SP_arguments: parse posicional
test/                       model_from_files, density_and_current_from_files, ...
samples/                    graphene, VSeTe ; test/ graphene, spin_graphene
```

### Objetos y firmas clave (verificadas)

```cpp
// hopping primitivo: (R[3], valor complejo, [i,j] base-0)
hopping_t = tuple<array<int,3>, complex<double>, array<int,2>>;

tuple<int /*num_wann*/, vector<string> /*lineas hopping*/>
read_wannier_file(const string);                 // DESCARTA ndegen (wz_gpoints)

hopping_list create_hopping_list(tuple<int, vector<string>>);  // 1-based->0-based, filtra ~0

hopping_list wrap_in_supercell(cellID_t cellDim, hopping_list); // tiling + PBC mod-wrap; COLAPSA q
void         save_supercell_as_csr(cellID_t, const hopping_list&, string); // paso único, byte-estable
```

Los 3 únicos call-sites de esas funciones (tbmodel.hpp:69, model_from_files.cpp:38
y :84) **encadenan** `create_hopping_list(read_wannier_file(...))`; ninguno
inspecciona la tupla intermedia.

Operadores actuales en `tbmodel`: velocidad = `H_ij·(−i)Δr` (diagonal-posición,
correcta); spin por etiquetas `_s+_/_s-_` (solo colineal sin mezcla SOC); sin L.

CSR de salida: `dim nnz` / valores `re im` / índices de columna / punteros de fila.
Entrada: `LABEL_hr.dat`, `LABEL.uc` (vectores de red), `LABEL.xyz` (label+posición).

---

## 1. Invariantes de diseño (rigen todos los planes)

1. **Un único objeto primitivo**: operador en espacio real `O_ij(R)` como
   `hopping_list` puro. `H, V, S, L` son instancias de lo mismo.
2. **Motor de expansión intacto**: `hopping_list → wrap_in_supercell →
   save_supercell_as_csr` no se toca. Toda fuente nueva desemboca ahí. Es la pieza
   rápida que alimenta a KPM.
3. **Normalización única**: `H` y todo `O` se normalizan (ndegen/wsvec) y expanden
   **idénticamente**, por la misma función. Nunca se corrige dos veces.
4. **Sidecar físico separado**: metadatos (componente cartesiana, observable,
   unidades, procedencia, cotas espectrales) viven aparte del `hopping_list`;
   nunca entran al kernel numérico.
5. **Estrangulador**: cada plan añade al lado; no destruye lo que funciona. La
   whitelist de ficheros es mínima y la funcionalidad, única.

**Caveats físicos asumidos** (conocidos):
- Velocidad diagonal-posición desprecia la conexión de Berry; OK para transporte
  DC, insuficiente para óptica/Berry. Modo conmutador completo (`r.dat`/`tb.dat`)
  como hito posterior.
- Corriente de spin `½{v,S}` es la *convencional*; con SOC no se conserva (la
  conservada de Shi et al. añade torque-dipolo). Estándar en spin-Hall KPM.
- Ruta proyector: `.amn` no ortonormal en general → Löwdin para cuantitativo.
- Base de Wannier ortonormal ⇒ `S = I` ⇒ este pipeline ejercita el camino
  **ortogonal** de KPM, **no** el `MetricPolicy` (`S≠I`) de LinQT Fase 4. Ese
  golden requiere una fuente no-ortogonal aparte (LCAO/NAO).

**Determinismo para goldens (byte-equality)**: orden fijo de malla-k;
normalización fija de FT (`1/N_k`, signo `e^{−ik·R}`); ordenación canónica de `R`
y de hoppings; umbral de truncado explícito; cuidado con el orden de sumación en
la FT.

---

## 2. Grafo de dependencias y fases

```text
Fase 0 (sin W90):   P1 ndegen ──► P2 ingesta externa/TB ──► P3 J derivado ──► P4 proyecto/seedname
                       │
Fase A (W90 H):        ├──► P5 wsvec
                       └──► P6 cotas (a,b) + descriptor
Fase B (SOC):       P7 gauge + .spn → spin exacto   (requiere P1, P6 para .eig)
Fase C:             P8 .amn + L locales → momento angular   (requiere P7: motor de gauge/FT)
```

P1 es precondición de todo (goldens fiables). P2–P4 cierran la Fase 0. P5–P6 la
Fase A. P7 la Fase B. P8 la Fase C. El motor de gauge+FT de P7 lo reutiliza P8.

---

## Plan 1 — Normalización `ndegen` correcta

**Funcionalidad única**: el parser captura `ndegen` (hoy lo descarta) y
`create_hopping_list` divide cada hopping por `ndegen[R]`. Aplica a cualquier `O`
que entre por esa ruta. En `ndegen ≡ 1` es **no-op**.

**Whitelist**
```text
include/wannier_parser.hpp   src/wannier_parser.cpp
include/hopping_list.hpp      src/hopping_list.cpp
test/ndegen_normalization.cpp  test/CMakeLists.txt
```
No se tocan `tbmodel.hpp` ni los tests existentes (solo encadenan): cambiar las
dos firmas en tándem los deja intactos.

**Cambio de firmas (en tándem)**
```cpp
tuple<int, vector<int> /*ndegen*/, vector<string>> read_wannier_file(const string);
hopping_list create_hopping_list(tuple<int, vector<int>, vector<string>>);
```

**Mapeo `ndegen[R]`**: `ndegen[k]` es el peso del k-ésimo punto WS. En `hr.dat` las
líneas van en bloques de `num_wann²` por punto WS, en el mismo orden que `ndegen`
(verificado en VSeTe: 484 = 22² líneas/bloque). Entonces
`indice_bloque = ordinal_linea_cruda / num_wann²`, `value /= ndegen[indice_bloque]`.
Contar sobre **líneas crudas**, no sobre hoppings conservados (el filtro de ceros
no debe desincronizar el contador).

**Corrección de bug latente (WS block)**: el parser lee `nrpts/15 + 1` líneas de
degeneración cuando lo correcto es `ceil(nrpts/15)`. Difieren cuando `nrpts` es
múltiplo de 15 → consume la 1ª línea de hopping como ndegen y desincroniza todo.
Enmascarado por los datos actuales (grafeno `nrpts=1`, VSeTe `nrpts=1037`). Fix
robusto: leer enteros token a token hasta tener exactamente `nrpts`, sin contar
líneas; asertar `ndegen.size()==nrpts`.

**Tests**: (1) regresión no-op grafeno byte-idéntica; (2) `hr.dat` sintético con
un bloque `ndegen=2` → valor a la mitad; (3) hopping ~0 → contador no se
desincroniza; (4) `nrpts=15` → no se come la 1ª línea de hopping.

**No-destrucción (verificada)**: grafeno/spin_graphene `ndegen≡1` → byte-idéntico,
tests pasan. VSeTe (sample sin test) corrige su salida.

**Desbloquea**: `H` correcto desde cualquier `hr.dat` real; precondición de todo
golden posterior.

---

## Plan 2 — Ingesta de operadores externos (`hr.dat`) / TB puro

**Funcionalidad única**: cargar cualquier `O_ij(R)` en formato `hr.dat` como
`hopping_list` y expandirlo por el motor existente, sin generarlo internamente.
Habilita modelos TB a mano (Haldane, Kane–Mele, BHZ), operadores externos y
debugging. TB puro = `hr.dat` con `ndegen ≡ 1`; un solo parser, no un formato
nuevo.

**Whitelist**
```text
include/tbmodel.hpp        # readOperatorModel(file) -> reusa create_hopping_list(read_wannier_file(file))
src/main.cpp               # rama para ingerir operador externo en vez de generarlo
include/w2sp_arguments.hpp # flag p.ej. --op-file NAME PATH ; o modo --raw
test/external_operator.cpp test/CMakeLists.txt
```

**Detalles**
- `readOperatorModel` reutiliza la cadena `create_hopping_list(read_wannier_file)`
  (ya con ndegen del P1; en raw todos a 1).
- Para `H` puro solo se necesita `hr.dat` (sin `.uc`/`.xyz`); velocidad sí
  necesita posiciones. La CLI debe permitir expandir un operador externo aunque
  falten `.uc`/`.xyz` si no se piden operadores derivados de posición.
- El motor de expansión no distingue procedencia: el operador externo va por el
  mismo `save_supercell_as_csr`.

**Tests / goldens analíticos**: grafeno NN (ya existe), Haldane, Kane–Mele, BHZ.
Validación física: reconstruir `H(k)=Σ_R e^{ik·R} H(R)` para unos k y comparar con
las bandas analíticas; o diagonalizar un supercell pequeño denso y comparar.

**Desbloquea**: primeros goldens con respuesta cerrada conocida; entrada de
operadores arbitrarios para los planes siguientes.

---

## Plan 3 — Corriente de spin derivada `J = ½(V·S + S·V)`

**Funcionalidad única**: salida derivada opcional que forma el anticonmutador
**post-expansión** como producto sparse de Eigen, sobre los supercell ya envueltos.

**Whitelist**
```text
include/operator_algebra.hpp  # anticommutator(SparseMatrix_t, SparseMatrix_t)
src/hopping_list.cpp / include # factorizar "escribir SparseMatrix_t a CSR" (hoy embebido en save_hopping_list_as_csr)
src/main.cpp                   # rama de salida J derivada
include/w2sp_arguments.hpp     # solicitar J (p.ej. JS_X_Z)
test/spin_current_anticomm.cpp test/CMakeLists.txt
```

**Detalles**
- Operar sobre `SparseMatrix_t` ensamblada (post-`wrap`), no sobre `hopping_list`:
  ensamblar `V` y `S`, `J = 0.5*(V*S + S*V)`, exportar.
- `V_sc` ya lleva el `Δr = r_J + A·q − r_I` correcto tras `wrap_in_supercell`;
  por eso el anticonmutador se forma después de expandir, no como convolución en
  `R` primitivo (que reintroduciría el bookkeeping de imágenes a mano).
- Fill-in acotado (alcance `r₁+r₂`); `V,S` hermíticos ⇒ `J` hermítico (chequeo).
- Factorizar el bloque de escritura CSR que hoy vive dentro de
  `save_hopping_list_as_csr` para reutilizarlo con una `SparseMatrix_t` directa.

**Granularidad**: wannier2sparse emite los **primitivos** `H,V,S,L`; `J` es
derivado opcional. LinQT puede hacerlo **lazy**: `J|ψ⟩ = ½(V(S|ψ⟩)+S(V|ψ⟩))` sin
materializar.

**Tests / cross-check**: para el caso colineal, `J` por matmul debe coincidir con
la `createHoppingSpinCurrents_list` por etiquetas existente → validación cruzada
fuerte. Más hermiticidad y un valor conocido.

**Decisión pendiente**: en Fase B, ¿se retira la ruta de spin por etiquetas o se
mantiene junto a la exacta?

**Desbloquea**: goldens de corriente de spin self-contained.

---

## Plan 4 — Proyecto/seedname (descubrimiento por carpeta)

**Funcionalidad única**: localizar ficheros por convención dentro de un directorio
(seedname), en lugar de flags sueltos por fichero. La CLI posicional actual sigue
funcionando (no destructivo).

**Whitelist**
```text
include/w2sp_arguments.hpp   # --project DIR / --seed NAME, resolución de rutas
src/main.cpp                 # usar rutas resueltas
test/project_resolution.cpp  test/CMakeLists.txt
```

**Detalles**
- Abstracción proyecto: dado `DIR` + `seedname`, resolver `seedname_hr.dat`,
  `seedname.uc`, `seedname.xyz`, y más adelante `.spn/.amn/u.mat/u_dis.mat/.eig`.
- Capa por encima del parsing posicional; mantener compatibilidad con la invocación
  actual `LABEL N1 N2 N3 [OP...]`.

**Tests**: apuntar a `samples/` y resolver los ficheros de VSeTe.

**Desbloquea**: "darle acceso a una carpeta de W90" — requisito explícito; base
para los parsers W90 de las fases B/C.

---

## Plan 5 — `wsvec` / `use_ws_distance`

**Funcionalidad única**: corrección de distancia mínima `R_eff = R + T` con reparto
de peso, aplicada **uniformemente** a `H` y a todo `O`.

**Whitelist**
```text
include/wannier_parser.hpp  src/wannier_parser.cpp   # read_wsvec
src/hopping_list.cpp                                  # aplicar en/ tras create_hopping_list
include/tbmodel.hpp                                   # cablear el fichero opcional
test/wsvec_split.cpp        test/CMakeLists.txt
```

**Detalles**
- Parsear `seedname_wsvec.dat`: por cada `(R,i,j)` una lista de vectores `T`
  (cuenta `nT` y `nT` líneas `T1 T2 T3`). **Verificar el formato exacto contra la
  fuente/doc de W90** antes de implementar.
- `O_eff(R+T) = O_raw(R,i,j) / (ndegen[R] · nT)`. Combinar con la división ndegen
  del P1 en la misma normalización.
- Para el golden, fijar `use_ws_distance` igual para todo en la corrida W90 (o
  ninguno).

**Tests**: `wsvec` sintético que parte un hopping en 2 `T` → mitad de peso y
desplazamiento correcto.

**Desbloquea**: `H` correcto para modelos modernos de W90 (uso común de
`use_ws_distance`).

---

## Plan 6 — Cotas espectrales `(a,b)` + descriptor (sidecar)

**Funcionalidad única**: emitir cotas espectrales y metadatos del operador en un
**descriptor separado** del kernel numérico.

**Whitelist**
```text
include/descriptor.hpp        # struct sidecar (componente, observable, unidades, procedencia, (a,b))
include/wannier_parser.hpp  src/wannier_parser.cpp   # read_eig (opcional)
src/main.cpp                  # emitir sidecar junto al CSR
test/spectral_bounds.cpp      test/CMakeLists.txt
```

**Detalles**
- `(E_min,E_max)` desde `.eig` (min/max sobre la malla-k) o unas pocas iteraciones
  de Lanczos/potencia sobre la `H` ensamblada.
- El descriptor no entra nunca en la recursión de Chebyshev (invariante 4).
- Encaja con el principio "bounds en el descriptor del operador" de LinQT y da
  procedencia al reescalado KPM.

**Tests**: cotas de un modelo analítico coinciden con la forma cerrada.

**Desbloquea**: procedencia espectral para el reescalado de Chebyshev en LinQT.

---

## Plan 7 — Motor de gauge + `.spn` → spin exacto (SOC)

**Funcionalidad única**: operador de spin exacto vía transformación de gauge:
`V(k)=U_dis(k)·U(k)`, `S_W(k)=V†S_B(k)V`, FT a `S_W(R)`. Sustituye la ruta de
etiquetas para modelos con mezcla SOC.

**Whitelist**
```text
include/gauge.hpp  src/gauge.cpp        # V=U_dis·U ; O_W(k)=V† O_B(k) V ; FT a O_W(R)
include/wannier_parser.hpp  src/wannier_parser.cpp   # read_umat, read_udis, read_spn (formateado), read_eig (de P6)
src/main.cpp  include/w2sp_arguments.hpp # ruta de spin exacto
test/gauge_spin.cpp  test/CMakeLists.txt
```

**Detalles**
- `.spn` **formateado** (`spn_formatted=.true.` en `pw2wannier90`), texto:
  `nbnd nkpts`, y por k el triángulo superior empaquetado de las 3 componentes
  (S hermítica: `S^a(n,m)=conj(S^a(m,n))`). **Verificar el orden exacto del
  empaquetado contra el fuente de `pw2wannier90`, no de memoria.**
- `V(k)=U_dis(k)·U(k)` (dims `num_bands×num_wann`); `S_W(k)=V†S_B(k)V`
  (`num_wann×num_wann`); `S_W(R)=(1/N_k)Σ_k e^{−ik·R}S_W(k)` sobre la **misma malla
  WS** que `hr.dat`. Salida como ficheros `O_ij(R)` → mismo motor.
- **Bookkeeping del disentanglement** (donde se rompen estos pipelines): las bandas
  de `.spn`, `.eig` y las filas de `u_dis.mat` deben referirse al mismo set,
  respetando `exclude_bands` y la ventana exterior. `S_W=V†S_B V` solo vale si
  `S_B(k)` está en el set sobre el que actúa `U_dis`. Asertar dimensiones por k.
- `u.mat`/`u_dis.mat` en **column-major Fortran**; malla k = `mp_grid` del `.win`.

**Tests**: modelo SOC pequeño de referencia; hermiticidad, reglas de suma;
comparar con WannierBerri o referencia si está disponible. Difícil de hacer
golden totalmente self-contained sin datos externos — anotar.

**Desbloquea**: `S` exacto y corriente de spin exacta (P3 con esta `S`) para SOC.

---

## Plan 8 — `.amn` + L locales → momento angular orbital

**Funcionalidad única**: momento angular orbital vía ruta proyector:
`C(k)=A†(k)V(k)`, `L_W(k)=C†L_local C`, FT a `L_W(R)`. Reutiliza el motor de
gauge/FT del P7.

**Whitelist**
```text
include/wannier_parser.hpp  src/wannier_parser.cpp   # read_amn
include/local_operators.hpp # Lx,Ly,Lz en base de armónicos reales p/d ; metadatos de proyección del .win
src/gauge.cpp / include      # ruta proyector C=A†V ; O_W(k)=C† O_local C (reutiliza P7)
src/main.cpp  include/w2sp_arguments.hpp
test/orbital_L.cpp  test/CMakeLists.txt
```

**Detalles**
- `A_mα(k)=⟨ψ_mk|g_α⟩` de `.amn`; `C(k)=A†(k)V(k)` (V de P7).
- `L_local` en armónicos reales p/d (matrices explícitas por shell); `L_W(k)=C†
  L_local C`; FT a `L_W(R)`.
- **Ortonormalización de Löwdin** de `A` si se quieren texturas cuantitativas.
- Metadatos de proyección (átomos/shells) del `.win`.
- Variante: spin local sobre átomos seleccionados = misma ruta proyector con
  `σ/2` en lugar de `L`.

**Tests**: autovalores `L_z` conocidos para un shell p o d; hermiticidad.

**Desbloquea**: operador de momento angular; corriente orbital (análogo de P3 con
`L`).

---

## 3. Decisiones pendientes (transversales)

1. **ABI**: el cambio de firma de P1 rompe la ABI si hubiera consumidores externos
   que enlacen `wannierlib` con la firma vieja. Asumido que no los hay (standalone,
   sale de LinQT). Si los hubiera → overload viejo.
2. **`hr.dat` podado**: el mapeo posicional `ndegen[bloque]` asume bloques
   completos `num_wann²` (canónico W90, verificado). Un `hr.dat` "podado" (solo
   no-ceros) rompería el posicional → mapa explícito `R→ndegen` leyendo el primer
   registro de cada bloque.
3. **Ruta de spin por etiquetas**: ¿se retira en P7 o convive con la exacta?
4. **Golden no-ortogonal**: `S≠I` (LinQT Fase 4) no lo da Wannier; ¿entra en
   alcance con otra fuente (LCAO/NAO) o se deja fuera?
