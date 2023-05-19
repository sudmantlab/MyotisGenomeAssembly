# Weird Genes to look at more carefully

## mMyoLuc1

### Miniprot

#### 10 largest genes:

N |  ID  |  length (bp)
------------------------
1 | UPI000216BBFE_1_431 | 1988265
2 | UPI00187CAF54_1_391 | 1456176
3 | UPI00178FFD5A_1_267 | 1454697
4 | UPI000216ED58_1_313 | 1268159
5 | UPI0000190B2E_1_64  | 1220555
6 | UPI0002A8980D_65_396 | 1206817
7 | UPI0017F498AF_29_114 | 1163911
8 | UPI0007040DE3_1_286 | 1142951
9 | UPI001852353E_1_107 | 1110772
10 | UPI000788203C_6_136 | 105516


### TOGA

#### 10 largest genes:

N |  ID  |  length (bp)
------------------------
1 | reg_13779 | 1570919
2 | reg_17759 | 1182896
3 | reg_2329 | 1050615
4 | reg_13260 1028547
5 | reg_13021 | 1009753
6 | reg_3877 | 966302
7 | reg_10575 | 942923
8 | reg_8671 | 915414
9 | reg_6383 | 848763
10 | reg_14521 | 848044

Its interesting to note that for reg_13779, its one mRNA that spans 1.6Mb. 
By grouping the GFF by genes, then calculating the fold difference between 
the shortest and longest isoforms, I can see how much of a problem this is.

**Note: This analysis is all based on the raw length of the mRNA feature in the GFF.**

`tibble(count = toga.fold.hist$count, bin = toga.fold.hist$mids)`

|   | count | bin  |
|------------------|
|   | <int> | <dbl>|
|1  | 25093 | 1000 |
|2  |   13  | 3000 |
|3  |    4  | 5000 |
|4  |    1  | 7000 |
|5  |    0  | 9000 |
|6  |    1  | 11000|
|7  |    0  | 13000|
|8  |    0  | 15000|
|9  |    1  | 17000|
|10 |    0  | 19000|
|11 |    0  | 21000|
|12 |    1  | 23000|

`tibble(count = toga.fold.under2k.hist$count, min = toga.fold.under2k.hist$breaks[1:20], max = toga.fold.under2k.hist$breaks[2:21])`

| count|  min|  max|
|-----:|----:|----:|
| 24733|    0|  100|
|   147|  100|  200|
|    74|  200|  300|
|    42|  300|  400|
|    20|  400|  500|
|    22|  500|  600|
|     9|  600|  700|
|    10|  700|  800|
|     8|  800|  900|
|     4|  900| 1000|
|     4| 1000| 1100|
|     3| 1100| 1200|
|     4| 1200| 1300|
|     3| 1300| 1400|
|     2| 1400| 1500|
|     1| 1500| 1600|
|     3| 1600| 1700|
|     1| 1700| 1800|
|     2| 1800| 1900|
|     1| 1900| 2000|

For comparison (using NCBI-based myoMyo annotations):

`tibble(count = myomyo.fold.hist$count, min = myomyo.fold.hist$breaks[1:20], max = myomyo.fold.hist$breaks[2:21]) %>% knitr::kable()`

| count| min| max|
|-----:|---:|---:|
| 19925|   0|  20|
|    23|  20|  40|
|     3|  40|  60|
|     1|  60|  80|
|     0|  80| 100|
|     0| 100| 120|
|     0| 120| 140|
|     0| 140| 160|
|     0| 160| 180|
|     0| 180| 200|
|     0| 200| 220|
|     0| 220| 240|
|     0| 240| 260|
|     0| 260| 280|
|     0| 280| 300|
|     0| 300| 320|
|     0| 320| 340|
|     0| 340| 360|
|     0| 360| 380|
|     1| 380| 400|

*** Very interesting to note that really, no myoMyo gene has a fold difference of over 20 - 
this is after the final product however. ***


Miniprot: 

`tibble(count = miniprot.fold.hist$count, min = miniprot.fold.hist$breaks[1:17], max = miniprot.fold.hist$breaks[2:18]) %>% knitr::kable()`

| count|  min|  max|
|-----:|----:|----:|
| 25978|    0|  500|
|   293|  500| 1000|
|    78| 1000| 1500|
|    36| 1500| 2000|
|    22| 2000| 2500|
|    12| 2500| 3000|
|    10| 3000| 3500|
|     5| 3500| 4000|
|     2| 4000| 4500|
|     2| 4500| 5000|
|     2| 5000| 5500|
|     0| 5500| 6000|
|     2| 6000| 6500|
|     0| 6500| 7000|
|     0| 7000| 7500|
|     1| 7500| 8000|
|     1| 8000| 8500|

Under 500:

| count| min| max|
|-----:|---:|---:|
| 22729|   0|  50|
|  1437|  50| 100|
|   589| 100| 150|
|   361| 150| 200|
|   285| 200| 250|
|   190| 250| 300|
|   117| 300| 350|
|   119| 350| 400|
|    80| 400| 450|
|    71| 450| 500|


Liftoff (myoMyo GenBank):

`tibble(count = liftoff.fold.hist$count, min = liftoff.fold.hist$breaks[1:20], max = liftoff.fold.hist$breaks[2:21]) %>% knitr::kable()`

| count|  min|  max|
|-----:|----:|----:|
| 20904|    0|  500|
|     0|  500| 1000|
|     0| 1000| 1500|
|     0| 1500| 2000|
|     0| 2000| 2500|
|     0| 2500| 3000|
|     0| 3000| 3500|
|     0| 3500| 4000|
|     0| 4000| 4500|
|     0| 4500| 5000|
|     0| 5000| 5500|
|     0| 5500| 6000|
|     0| 6000| 6500|
|     0| 6500| 7000|
|     1| 7000| 7500|

And those under 500:

| count| min| max|
|-----:|---:|---:|
| 20484|   0|  10|
|   266|  10|  20|
|    72|  20|  30|
|    31|  30|  40|
|    18|  40|  50|
|    10|  50|  60|
|     4|  60|  70|
|     2|  70|  80|
|     4|  80|  90|
|     2|  90| 100|
|     2| 100| 110|
|     1| 110| 120|
|     3| 120| 130|
|     1| 130| 140|
|     0| 140| 150|
|     0| 150| 160|
|     0| 160| 170|
|     0| 170| 180|
|     1| 180| 190|
|     2| 190| 200|
|     1| 200| 210|

*** When you compare the genes with the top fold difference in mRNA sizes between myoMyo's 
original annotation and its liftoff, they are mutually exclusive in the top 50. ***

ATXN1 and DMD are in the top 50 in the original myoMyo (as you'd expect), but not in the liftoff.

ATXN1 does have a similar fold difference between isoforms (~19) in both, but DMD has a 
fold difference of ~1 for the liftoff, and of ~32 for the original. 


