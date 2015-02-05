# Searching GEO
Brian High  
2/3/2015  

## Setting up some options

Let's first turn on the cache for increased performance and improved styling.


```r
# Set some global knitr options
library("knitr")
opts_chunk$set(tidy=TRUE, tidy.opts=list(blank=FALSE, width.cutoff=60), 
               cache=FALSE, messages=FALSE)
```

Load the `pandoc` package so we can make nicer table listings with `pandoc.table`.


```r
suppressMessages(library(pander))
```

## Prepare for HW2

We will need to query the GEOmetadb database. Today we will explore this 
database and practice various ways to query it.

## Load the `GEOmetadb` package

First we load the `GEOmetadb` library.


```r
suppressMessages(library(GEOmetadb))
```

Let's also view the available methods.


```r
ls("package:GEOmetadb")
```

```
## [1] "columnDescriptions" "geoConvert"         "getBiocPlatformMap"
## [4] "getSQLiteFile"
```

## Download the GEO database

We should have already downloaded this database when viewing the lecture slides.


```r
## This will download the entire database, so can be slow
if (!file.exists("GEOmetadb.sqlite")) {
    # Download database only if it's not done already
    getSQLiteFile()
}
```

## List tables with `SQL`

In `SQL`, you can query the database structure with ordinary `SQL` commands.


```r
geo_con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
dbGetQuery(geo_con, "SELECT name FROM sqlite_master WHERE type='table';")
```

```
##                 name
## 1                gse
## 2                gpl
## 3                gsm
## 4            gse_gsm
## 5            gse_gpl
## 6                gds
## 7         gds_subset
## 8            sMatrix
## 9  geodb_column_desc
## 10        geoConvert
## 11          metaInfo
```

## List `gse` fields with `SQL`

The `PRAGMA` command is a standard `SQLite` command.


```r
dbGetQuery(geo_con, "PRAGMA table_info(gse);")
```

```
##    cid                 name    type notnull dflt_value pk
## 1    0                   ID    REAL       0       <NA>  0
## 2    1                title    TEXT       0       <NA>  0
## 3    2                  gse    TEXT       0       <NA>  0
## 4    3               status    TEXT       0       <NA>  0
## 5    4      submission_date    TEXT       0       <NA>  0
## 6    5     last_update_date    TEXT       0       <NA>  0
## 7    6            pubmed_id INTEGER       0       <NA>  0
## 8    7              summary    TEXT       0       <NA>  0
## 9    8                 type    TEXT       0       <NA>  0
## 10   9          contributor    TEXT       0       <NA>  0
## 11  10             web_link    TEXT       0       <NA>  0
## 12  11       overall_design    TEXT       0       <NA>  0
## 13  12              repeats    TEXT       0       <NA>  0
## 14  13  repeats_sample_list    TEXT       0       <NA>  0
## 15  14             variable    TEXT       0       <NA>  0
## 16  15 variable_description    TEXT       0       <NA>  0
## 17  16              contact    TEXT       0       <NA>  0
## 18  17   supplementary_file    TEXT       0       <NA>  0
```

## List tables with `dbListTables`

Instead of using `SQL` commands, we can list tables and fields with functions
from the `GEOmetadb` package.


```r
geo_con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
dbListTables(geo_con)
```

```
##  [1] "gds"               "gds_subset"        "geoConvert"       
##  [4] "geodb_column_desc" "gpl"               "gse"              
##  [7] "gse_gpl"           "gse_gsm"           "gsm"              
## [10] "metaInfo"          "sMatrix"
```


```r
dbListFields(geo_con, "gse")
```

```
##  [1] "ID"                   "title"                "gse"                 
##  [4] "status"               "submission_date"      "last_update_date"    
##  [7] "pubmed_id"            "summary"              "type"                
## [10] "contributor"          "web_link"             "overall_design"      
## [13] "repeats"              "repeats_sample_list"  "variable"            
## [16] "variable_description" "contact"              "supplementary_file"
```

## Explore `gse`


```r
columnDescriptions()[1:5, ]
```

```
##   TableName        FieldName
## 1       gse            title
## 2       gse              gse
## 3       gse           status
## 4       gse  submission_date
## 5       gse last_update_date
##                                                Description
## 1                 unique name describing the overall study
## 2 unique accession number approved and issued by GEO, NCBI
## 3                                  date released to public
## 4                                           date submitted
## 5                                        date last updated
```

## Load library `data.table`

This will provide us with some practice querying with data.table.


```r
suppressMessages(library(data.table))
```

## Explore `gse` with `data.table`


```r
cd <- as.data.table(columnDescriptions())
cd[TableName == "gse", FieldName]
```

```
##  [1] "title"                "gse"                  "status"              
##  [4] "submission_date"      "last_update_date"     "pubmed_id"           
##  [7] "summary"              "type"                 "contributor"         
## [10] "contact"              "web_link"             "overall_design"      
## [13] "repeats"              "repeats_sample_list"  "variable"            
## [16] "variable_description" "supplementary_file"
```

## List `gse` columns with `pandoc.table`


```r
gsefields <- as.data.frame(cd[TableName == "gse" & FieldName %in% 
    c("gse", "title", "pubmed_id", "summary", "contact")])
pandoc.table(gsefields, style = "grid")
```

```
## 
## 
## +-------------+-------------+--------------------------------+
## |  TableName  |  FieldName  |          Description           |
## +=============+=============+================================+
## |     gse     |    title    |   unique name describing the   |
## |             |             |         overall study          |
## +-------------+-------------+--------------------------------+
## |     gse     |     gse     |    unique accession number     |
## |             |             |  approved and issued by GEO,   |
## |             |             |              NCBI              |
## +-------------+-------------+--------------------------------+
## |     gse     |  pubmed_id  |  [Values separated by ';tab',  |
## |             |             | if more than one] NCBI PubMed  |
## |             |             |       identifier (PMID)        |
## +-------------+-------------+--------------------------------+
## |     gse     |   summary   | a description of the goals and |
## |             |             |    objectives of this study    |
## +-------------+-------------+--------------------------------+
## |     gse     |   contact   |  contact information for this  |
## |             |             |             study              |
## +-------------+-------------+--------------------------------+
```

## Explore `gpl`


```r
cd[TableName == "gpl", FieldName]
```

```
##  [1] "ID"                   "title"                "gpl"                 
##  [4] "status"               "submission_date"      "last_update_date"    
##  [7] "technology"           "distribution"         "organism"            
## [10] "manufacturer"         "manufacture_protocol" "coating"             
## [13] "catalog_number"       "support"              "description"         
## [16] "web_link"             "contact"              "data_row_count"      
## [19] "supplementary_file"   "bioc_package"
```

## Explore columns in `gpl`


```r
gplfields <- as.data.frame(cd[TableName == "gpl" & FieldName %in% 
    c("gpl", "organism", "manufacturer")])
pandoc.table(gplfields, style = "grid")
```

```
## 
## 
## +-------------+--------------+-------------------------------+
## |  TableName  |  FieldName   |          Description          |
## +=============+==============+===============================+
## |     gpl     |     gpl      | unique GEO Platfrom accession |
## |             |              | number approved and issued by |
## |             |              |           GEO, NCBI           |
## +-------------+--------------+-------------------------------+
## |     gpl     |   organism   | [Values separated by ';tab',  |
## |             |              | if more than one] Organism(s) |
## +-------------+--------------+-------------------------------+
## |     gpl     | manufacturer | name of the company, facility |
## |             |              | or laboratory where the array |
## |             |              | was manufactured or produced  |
## +-------------+--------------+-------------------------------+
```

## Explore `gse_gpl`


```r
cd[TableName == "gse_gpl", FieldName]
```

```
## [1] "gse" "gpl"
```

## Explore columns in `gse_gpl`

Why are there only two fields in this table? What is this table for?


```r
gse_gplfields <- as.data.frame(cd[TableName == "gse_gpl"])
pandoc.table(gse_gplfields, style = "grid")
```

```
## 
## 
## +-------------+-------------+-------------------+
## |  TableName  |  FieldName  |    Description    |
## +=============+=============+===================+
## |   gse_gpl   |     gse     |  GEO Series name  |
## +-------------+-------------+-------------------+
## |   gse_gpl   |     gpl     | GEO Platform name |
## +-------------+-------------+-------------------+
```

## List "title" fields with `pandoc.table`

Why do many tables include a "title" field? Are the titles the same?


```r
gsefields <- as.data.frame(cd[FieldName == "title"])
pandoc.table(gsefields, style = "grid")
```

```
## 
## 
## +-------------+-------------+-----------------------------+
## |  TableName  |  FieldName  |         Description         |
## +=============+=============+=============================+
## |     gse     |    title    | unique name describing the  |
## |             |             |        overall study        |
## +-------------+-------------+-----------------------------+
## |     gpl     |    title    | unique name describing the  |
## |             |             |       Platform (GPL)        |
## +-------------+-------------+-----------------------------+
## |     gsm     |    title    | unique name describing this |
## |             |             |           Sample            |
## +-------------+-------------+-----------------------------+
## |     gds     |    title    |      title of this GDS      |
## +-------------+-------------+-----------------------------+
```

## List "contact" field structure

Let's look at some records in `gse`. What does a "contact" look like?


```r
query <- "SELECT contact FROM gse LIMIT 1;"
res <- dbGetQuery(geo_con, query)
strsplit(res$contact, "\t")
```

```
## [[1]]
##  [1] "Name: Michael Bittner;"                                                     
##  [2] "Email: mbittner@nhgri.nih.gov;"                                             
##  [3] "Phone: 301-496-7980;"                                                       
##  [4] "Fax: 301-402-3241;"                                                         
##  [5] "Department: Cancer Genetics Branch;"                                        
##  [6] "Institute: NHGRI, NIH;"                                                     
##  [7] "Address:  ;"                                                                
##  [8] "City: Bethesda;"                                                            
##  [9] "State: MD;"                                                                 
## [10] "Zip/postal_code: 20892;"                                                    
## [11] "Country: USA;"                                                              
## [12] "Web_link: http://www.nhgri.nih.gov/Intramural_research/People/bittnerm.html"
```

## Find manufacturer data

Query the manufacturers with a `SQL` command, listed with `data.table`...


```r
manu <- data.table(dbGetQuery(geo_con, 
    "SELECT DISTINCT manufacturer FROM gpl ORDER BY manufacturer ASC;"))
manu[,list(length(manufacturer)), by=manufacturer]
```

```
##                                         manufacturer V1
##    1:                                             NA  1
##    2:                                              -  1
##    3:                                              .  1
##    4:                                            454  1
##    5:                              454 Life Sciences  1
##   ---                                                  
## 2075: washington university microarray core facility  1
## 2076:         www.MYcroarray.com, Ann Arbor, MI, USA  1
## 2077:                                www.agilent.com  1
## 2078:                           www.chem.agilent.com  1
## 2079:                            www.combimatrix.com  1
```

## Our `SQL` command

We just wanted a list of manufacturers so the `SQL` query is:

```
SELECT DISTINCT manufacturer FROM gpl 
ORDER BY manufacturer ASC;
```

However, since we also grouped `by=manufacturer` in our `data.table`, we could 
have simply used the `SQL` query:

```
SELECT manufacturer FROM gpl;
```
Let's try that...

## Find manufacturer data

Query the manufacturers with a simpler `SQL` command ... grouping with `by` and 
ordering with `setkey` in `data.table`...


```r
manu <- data.table(dbGetQuery(geo_con, 
            "SELECT manufacturer FROM gpl;"))
setkey(manu, manufacturer)
manu[,list(length(manufacturer)), by=manufacturer]
```

```
##                                         manufacturer V1
##    1:                                             NA  1
##    2:                                              -  1
##    3:                                              .  1
##    4:                                            454  1
##    5:                              454 Life Sciences  1
##   ---                                                  
## 2075: washington university microarray core facility  1
## 2076:         www.MYcroarray.com, Ann Arbor, MI, USA  1
## 2077:                                www.agilent.com  1
## 2078:                           www.chem.agilent.com  1
## 2079:                            www.combimatrix.com  1
```


## Finding data with a `join`

To get supplementary file names ending with `CEL.gz` (case-insensitive) from 
only manufacturer Affymetrix, we need to `join` the `gsm` and `gpl` tables. 

```
SELECT 
        gpl.bioc_package, 
        gsm.title, 
        gsm.series_id, 
        gsm.gpl, 
        gsm.supplementary_file 
    FROM gsm 
    JOIN gpl ON gsm.gpl=gpl.gpl 
    WHERE gpl.manufacturer='Affymetrix' 
        AND gsm.supplementary_file like '%CEL.gz';
```

## Now let's run that query


```r
query<-"SELECT 
            gpl.bioc_package, 
            gsm.title, 
            gsm.series_id, 
            gsm.gpl, 
            gsm.supplementary_file 
        FROM gsm 
        JOIN gpl ON gsm.gpl=gpl.gpl 
        WHERE gpl.manufacturer='Affymetrix' 
            AND gsm.supplementary_file like '%CEL.gz';"
res <- dbGetQuery(geo_con, query)
head(res, 3)
```

```
##   bioc_package      title series_id   gpl
## 1       hu6800 BM_CD34-1a    GSE500 GPL80
## 2       hu6800 BM_CD34-1b    GSE500 GPL80
## 3       hu6800  BM_CD34-2    GSE500 GPL80
##                                                         supplementary_file
## 1 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnn/GSM575/suppl/GSM575.cel.gz
## 2 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnn/GSM576/suppl/GSM576.cel.gz
## 3 ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnn/GSM577/suppl/GSM577.cel.gz
```

## Why did we need a `join`?

The 
[GEOmetadb database](http://gbnci.abcc.ncifcrf.gov/geo/geo_help.php), 
is a [relational database](http://en.wikipedia.org/wiki/Relational_database). 

There are several tables which can be linked on common fields. 

Since each table 
contains data for only one type of record, tables must be linked to search for 
fields pertaining to the various types of records. 

We join on the common fields, 
called [keys](http://en.wikipedia.org/wiki/Relational_database#Primary_key).

## Table Relationships of `GEOmetadb`

![Table Relationships](http://gbnci.abcc.ncifcrf.gov/geo/images/GEOmetadb_diagram.png)

Source: [Help: GEOmetadb Application, Meltzerlab/GB/CCR/NCI/NIH &copy;2008](http://gbnci.abcc.ncifcrf.gov/geo/geo_help.php)

## Keys of `GEOmetadb`

```
+------------+-------+------------------------------------------------+
| Table      | Key   | Links to Table.Key                             |
+============+=======+================================================+
| gse        | gse   | gse_gpl.gse, gse_gsm.gse, gds.gse, sMatrix.gse |
+------------+-------+------------------------------------------------+
| gpl        | gpl   | gds.gpl, gse_gpl.gpl, sMatrix.gpl, gsm.gpl     |
+------------+-------+------------------------------------------------+
| gsm        | gsm   | gse_gsm.gsm                                    |
| gsm        | gpl   | gds.gpl, gse_gpl.gpl, sMatrix.gpl, gpl.gpl     |
+------------+-------+------------------------------------------------+
| gds        | gds   | gds_subset.gds                                 |
+------------+-------+------------------------------------------------+
| gds_subset | gds   | gds.gds                                        |
+------------+-------+------------------------------------------------+
| sMatrix    | gse   | gse_gpl.gse, gse_gsm.gse, gds.gse, gse.gse     |
| sMatrix    | gpl   | gds.gpl, gse_gpl.gpl, gpl.gpl, gsm.gpl         |
+------------+-------+------------------------------------------------+
| gse_gpl    | gse   | gse_gpl.gse, gse_gsm.gse, gds.gse, sMatrix.gse |
| gse_gpl    | gpl   | gds.gpl, gse_gpl.gpl, gpl.gpl, sMatrix.gpl     |
+------------+-------+------------------------------------------------+
| gse_gsm    | gse   | gse_gpl.gse, gse.gse, gds.gse, sMatrix.gse     |
| gse_gsm    | gsm   | gsm.gsm                                        |
+------------+-------+------------------------------------------------+
```

Source: [Help: GEOmetadb Application, Meltzerlab/GB/CCR/NCI/NIH &copy;2008](http://gbnci.abcc.ncifcrf.gov/geo/geo_help.php)

## A three-table `join`

To get raw data, we need to `join` three tables with two `join` clauses...


```r
query<-"SELECT gsm.gsm, gsm.supplementary_file 
        FROM (gse JOIN gse_gsm ON gse.gse=gse_gsm.gse) j 
        JOIN gsm ON j.gsm=gsm.gsm 
        WHERE gse.pubmed_id='21743478' 
        LIMIT 2;"
res <- as.data.table(dbGetQuery(geo_con, query))
res[,strsplit(gsm.supplementary_file, ';\t'),by=gsm.gsm]
```

```
##      gsm.gsm
## 1: GSM733816
## 2: GSM733816
## 3: GSM733817
## 4: GSM733817
##                                                                                   V1
## 1: ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733816/suppl/GSM733816.CEL.gz
## 2: ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733816/suppl/GSM733816.chp.gz
## 3: ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733817/suppl/GSM733817.CEL.gz
## 4: ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM733nnn/GSM733817/suppl/GSM733817.chp.gz
```
