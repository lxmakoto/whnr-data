# WhnR — datos de secuencia, alineamientos y análisis de selección

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23439194.svg)](https://doi.org/10.5281/zenodo.23439194)

Datos que acompañan al manuscrito sobre **WhnR** (`whnR`, antes
`cl33274978` en AgnostosDB), una fusión winged-helix / nucleotidiltransferasa
tipo pol-β conservada en el superfilo Patescibacteria.

Autor: Oskar Kruk · `oscar.kruk@autonoma.cat`

Cítese como: Kruk, O. *Datos de WhnR: secuencias codificantes, alineamientos y análisis de selección*. Zenodo. [doi:10.5281/zenodo.23439194](https://doi.org/10.5281/zenodo.23439194)

---

## `secuencias/` — secuencias codificantes

Las CDS nucleotídicas de los homólogos, recuperadas siguiendo la cadena
accesión → UPI de UniParc → referencia cruzada EMBL/EMBLWGS → `protein_id` →
registro nucleotídico original.

| Fichero | Contenido |
|---|---|
| `familia_codones.fasta` | **716 CDS** alineadas en codones, 573 nt (191 codones). De 742 homólogos, 740 tenían CDS recuperable y 716 traducen exactamente a su fila del alineamiento |
| `cl33274978_codones.fasta` | 27 CDS extraídas directamente de los MAG por coordenadas, 549 nt |
| `cds_cache.json` | Caché de la recuperación: accesión → UPI → `protein_id` → secuencia |
| `pid_map.json` | Correspondencia accesión → `protein_id` |

**Verificación:** las 27 CDS de los MAG traducen exactamente al producto
almacenado, y ningún codón de ninguno de los dos alineamientos discrepa del
aminoácido de su columna.

## `alineamientos/`

| Fichero | Contenido |
|---|---|
| `familia_742.a3m` | Alineamiento de la familia: 742 homólogos, servidor MSA de ColabFold (UniRef30 + ColabFoldDB) |
| `control_KNTasa_1469.a3m` | **Control emparejado**: 1469 homólogos de la nucleotidiltransferasa de kanamicina (P05057), del mismo servidor, para comparar conservación del par catalítico con el mismo código |
| `cl33274978_proteina.fasta` | Alineamiento de proteína de los 27 miembros del clúster |
| `cl33274978_codones_dedup99.fasta` | 8 genotipos tras colapsar casi clones al 99 % (clado Doudnabacteria) |
| `subfamilia_codones.fasta` | 15 secuencias de un clado independiente (Moranbacteria / Parcubacteria) |

## `datamonkey/`

Entradas y salidas de los análisis de selección en
[Datamonkey](https://datamonkey.org).

**Entradas** (alineamiento de codones con el árbol Newick incrustado):

- `DATAMONKEY_dedup99.fasta` — 8 genotipos, clado Doudnabacteria
- `DATAMONKEY_subfamilia.fasta` — 15 secuencias, clado Moranbacteria
- `DATAMONKEY_completo.fasta` — los 27 MAG, **solo para comparación**: contiene
  casi clones y produce soporte inflado

**Salidas** (JSON crudo del servidor):

| Fichero | Método | Conjunto |
|---|---|---|
| `fel_results.json` | FEL | Doudnabacteria |
| `fubar_results.json` | FUBAR | Doudnabacteria |
| `meme_results.json` | MEME | Doudnabacteria |
| `fel_subfamilia.json` | FEL | Moranbacteria |
| `meme_subfamilia.json` | MEME | Moranbacteria |

## `tablas/`

| Fichero | Contenido |
|---|---|
| `seleccion_por_sitio.tsv` | Los 183 sitios con conservación ponderada, FEL, FUBAR y MEME del clado Doudnabacteria |
| `fel_subfamilia_por_sitio.tsv` | FEL del clado Moranbacteria, comparado con el anterior |
| `meme_subfamilia.tsv` | MEME de Moranbacteria con corrección FDR |
| `fel_por_sitio.tsv` | FEL de Doudnabacteria, tabla completa |
| `mapa_columnas.tsv` | Columna del alineamiento → número de residuo de la referencia |
| `ntase_activas.json` | 21 nucleotidiltransferasas pol-β activas con estructura resuelta: par catalítico y distancias entre carboxilatos |

## `arboles/`

Árboles de unión de vecinos en formato Newick. `_k2p` usa distancias corregidas
por Kimura 2 parámetros; los demás, distancias p sin corregir. **La topología es
idéntica entre ambos** (Robinson-Foulds = 0).

---

## Numeración

Los alineamientos de codones del clúster cubren los **residuos 1–183** de los
191 de la referencia (`GCA_001777825.1`, el representante de AgnostosDB), con
correspondencia **1:1**: columna *n* = residuo *n*. Los alineamientos de la
familia cubren los 191.

Residuos relevantes:

| Residuo | Conservación en 742 homólogos | Comentario |
|---|---|---|
| **Asp121** | 1,00 | Enterrado (SASA 0,066); purificadora en los dos clados |
| **Asp166** | 0,99 | Puente salino con Arg161 (2,87 Å) |
| **Arg161** | 0,99 | |
| Glu119 | 0,18 como D/E | Posición *i* del par catalítico: **no es carboxilato**. Glu/Pro/Phe en Doudnabacteria, Lys en Moranbacteria |

## Advertencias de uso

1. **Los sitios sinónimos están saturados a escala de superfilo.** La identidad
   en 3.ª posición del codón de `familia_codones.fasta` es de 0,306 frente a
   ≈0,27 al azar. Ese fichero **no sirve para dN/dS**; se incluye por su valor
   como recurso de secuencia. Usa los subconjuntos por clado.
2. El subconjunto de Doudnabacteria (0,347) está en el límite de la saturación.
   El de Moranbacteria (0,473) es el de mejor señal.
3. Los árboles son de unión de vecinos, no de máxima verosimilitud.

## `scripts/`

Los scripts que generan todo lo anterior, en orden de ejecución.

| Script | Qué hace |
|---|---|
| `54_codon_aln.py` | Extrae las CDS de los 27 MAG por coordenadas, verifica la traducción y construye el alineamiento de codones |
| `55_datamonkey_extras.py` | Mapa columna→residuo, deduplicación al 99 % y árboles |
| `56_msa_control.py` | Pide al servidor MSA de ColabFold los homólogos de la nucleotidiltransferasa activa |
| `57_control_emparejado.py` | Conservación de ambas familias con el mismo código y ponderación |
| `58_cds_familia.py` | Recupera las CDS de los 742 homólogos vía UniParc → EMBL → registro nucleotídico |
| `59_codones_familia.py` | Verifica y retrotraduce las 716 CDS al alineamiento de codones |
| `60_subconjunto_dnds.py` | Primer intento de subconjunto no saturado (enlace simple; encadena y falla) |
| `61_subconjunto_v2.py` | Subconjunto compacto exigiendo que todos los pares superen el umbral |
| `63_distribucion_sitio.py` | Distancias entre carboxilatos en 21 pol-β activas |

`60` se incluye aunque no funcione: documenta por qué el enlace simple no sirve
para este problema.

Dependencias: `numpy`, `pandas`, `biopython`, `pyhmmer`, `matplotlib`.

---

## Licencia

[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/)
(CC BY 4.0), igual que el depósito de Zenodo. Cubre tanto los datos como los
scripts: puedes usarlos, modificarlos y redistribuirlos, también con fines
comerciales, citando la fuente.

Las secuencias de partida proceden de AgnostosDB, UniParc, ENA/GenBank y el PDB,
cada una con sus propios términos de uso.
