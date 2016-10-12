# Data structures in R

## Basic data structures



|           |dimensions |Data.types       |
|:----------|:----------|:----------------|
|vector     |1          |1                |
|matrix     |2          |1                |
|data.frame |2          |n (1 per column) |
|list       |1          |n                |
|array      |n          |1                |

## Specialised data structures

Such specialised data structures, also called *classes* are

- **defined** for specialised domains like proteomics, genomics
  (micro-arrays, high-throughput sequencing, ...)

- **composed** of basic data structures like those summarised above (and
  others) that store bits that form the specialised *class*
  
- All parts are **contained/unified** into one *class* and are
  manipulated/processed consistently.
  
We will see several of these proteomics *ad hoc* classes (for
quantitative data, for PSMs identification, for raw mass spectrometry
data) and how to use them throughout the workshop.

