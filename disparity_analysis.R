########################################################################
## Title: Poverty and POU disparity anlaysis.
## Date: 2014-04-21
########################################################################

## NOTE (Michael): Download more indicator from the world bank for
##                 imputation and analysis.
##

## NOTE (Michael): Need to have a transformation table, so I can see
##                 the difference of imputation.

## NOTE (Michael): Next step is to examine the distributional
##                 difference between the country with problem and
##                 without.

library(car)
library(zoo)
library(data.table)
library(lattice)
library(FAOSTAT)
library(Amelia)
library(ggplot2)
library(lme4)
refreshPovData = FALSE

## Data preperation step
## ---------------------------------------------------------------------
## Load the data and remove future years
load(file = "final.RData")
fsiSubYear.df = fsi.df[fsi.df$Year %in% c(1995:2010), ]

## Read in the meta-data
meta.dt = data.table(read.csv(file = "meta.csv",
    stringsAsFactors = FALSE))

## Country name table
countryNamesTable =
    data.table(FAOcountryProfile[, c("FAOST_CODE", "FAO_TABLE_NAME")])
countryNamesTable[FAOST_CODE == 357,
                    FAO_TABLE_NAME := "Taiwan and China"]
countryNamesTable[FAOST_CODE == 351,
                    FAO_TABLE_NAME := "China Aggregate"]
countryNamesTable[FAOST_CODE == 107, FAO_TABLE_NAME := "Cote d'Ivoire"]
countryNamesTable[FAOST_CODE == 284, FAO_TABLE_NAME := "Aland Islands"]
countryNamesTable[FAOST_CODE == 279, FAO_TABLE_NAME := "Curacao"]
countryNamesTable[FAOST_CODE == 182, FAO_TABLE_NAME := "Reunion"]
countryNamesTable[FAOST_CODE == 282,
                    FAO_TABLE_NAME := "Saint Barthelemy"]
countryNamesTable[FAOST_CODE == 206,
                    FAO_TABLE_NAME := "former Sudan"]

## Download poverty and GINI data
if(refreshPovData){
    pov.lst =
        getWDItoSYB(indicator = meta.dt[!SOFI & SOURCE == "raw",
                        VARIABLE])
    pov.dt = data.table(merge(pov.lst$entity,
        FAOcountryProfile[, c("ISO2_WB_CODE", "FAOST_CODE")],
        all.x = TRUE))
    pov.dt[, `:=`(ISO2_WB_CODE = NULL, Country = NULL)]
    ## pov.df$ISO2_WB_CODE = NULL
    ## pov.df$Country = NULL
    save(pov.dt, file = "poverty.Rdata")
}

load("poverty.Rdata")
     
## Final unimputed data set
final.dt =
    data.table(fsiSubYear.df[, c("FAOST_CODE", "Year",
                                 meta.dt[SOFI & SOURCE == "raw",
                                         VARIABLE])])
final.dt = merge(final.dt, pov.dt, by = c("FAOST_CODE", "Year"))
setkeyv(final.dt, c("FAOST_CODE", "Year"))
final.dt[, Year := as.integer(Year)]

## Check the sparsity
sapply(final.dt, FUN = function(x) round(sum(is.na(x))/length(x) * 100))

## Check variable type
sapply(final.dt, FUN = function(x) range(x, na.rm = TRUE))

## Make correction to irrigation share of arable land.
final.dt[RL.AREA.EQIRR.HA.SHL >= 100, RL.AREA.EQIRR.HA.SHL := 100]

## Data transformation step:
## ---------------------------------------------------------------------

## Scale proportion data
scaleString = paste0("`:=`(",
    paste0(meta.dt[SCALE != 1, VARIABLE], " = ",
           meta.dt[SCALE != 1, VARIABLE], "/",
           meta.dt[SCALE != 1, SCALE], collapse = ", "), ")")
final.dt[, eval(parse(text = scaleString))]


## Some data transformation for normality and interpretation
invlogit = function(x){
    x = ifelse(x == 1, 1 - 1e-5,
        ifelse(x == 0, 0 + 1e-5, x))
    log(x/(1 - x))
}

transformString = paste0("`:=`(",
    paste0(meta.dt[!is.na(TRANSFORM) & SOURCE == "raw", VARIABLE], " = ",
           ifelse(meta.dt[!is.na(TRANSFORM) & SOURCE == "raw",
                          TRANSFORM] == "logit", "invlogit(", "log("),
           meta.dt[!is.na(TRANSFORM) & SOURCE == "raw", VARIABLE], ")",
           collapse = ", "), ")")
final.dt[, eval(parse(text = transformString))]

## Removing strange countries
final.dt = final.dt[FAOST_CODE != 164, ]
final.dt = final.dt[FAOST_CODE != 186, ]

plotString =
    paste0("~ ", paste0(colnames(final.dt)[-c(1, 2)], collapse = " + "))
scatterplotMatrix(eval(parse(text = plotString)), data = final.dt,
                  smoother = FALSE, use = "pairwise.complete.obs")


## Imputation Step
## ---------------------------------------------------------------------

## Create bounds for imputation

variableType = sapply(final.dt, typeof)
boundColumns = which(variableType == "double")
boundUpper = sapply(final.dt[, names(boundColumns), with = FALSE],
    max, na.rm = TRUE) +
         2 * sapply(final.dt[, names(boundColumns), with = FALSE], sd,
                    na.rm = TRUE)
boundLower = sapply(final.dt[, names(boundColumns), with = FALSE],
    max, na.rm = TRUE) -
         2 * sapply(final.dt[, names(boundColumns), with = FALSE], sd,
                    na.rm = TRUE)
boundsMatrix = matrix(c(boundColumns, boundLower, boundUpper),
    byrow = FALSE, nc = 3)
rownames(boundsMatrix) = names(boundColumns)

## lower bounds for untransformed variable
boundsMatrix[meta.dt[SOURCE == "raw" &
                     TRANSFORM == "" &
                     TYPE == "positive", VARIABLE],
             2] = 0


## linear interpolation for available data
myApprox = function(x){
    n = length(na.omit(x))
    if(n >= 2){
        tmp = na.approx(x, na.rm = FALSE)
    } else {
        tmp = x
    }
    tmp
}

for(i in colnames(final.dt)[-c(1:2)]){
    final.dt[, eval(parse(text = paste0(i, " := myApprox(", i, ")"))),
             by = "FAOST_CODE"]
}

## multiply impute the data set
imputed.am = amelia(final.dt, ts = "Year", cs = "FAOST_CODE",
    lag = c("ADESA", "POPULATION"), 
    p2s = 1, bounds = boundsMatrix, m = 10)

## Combine the multiple imputation by fitting loess
mi = Reduce(f = function(x, y) rbind(x, y), x = imputed.am$imputations)
combineString =
    paste0("list(", paste0(paste0(meta.dt[SOURCE == "raw", VARIABLE],
    "= predict(loess(", meta.dt[SOURCE == "raw", VARIABLE],
    " ~ Year, span = 0.85, loess.control(surface = 'direct')),",
    "newdata = data.frame(Year = 1995:2010))"),
    collapse = ", "), ")")
imputed.dt = mi[, eval(parse(text = combineString)), by = "FAOST_CODE"]
imputed.dt[, Year := rep(1995:2010, length(unique(FAOST_CODE)))]

## Plot the relationship of the result
plotString =
    paste0("~ ", paste0(colnames(imputed.dt)[colnames(imputed.dt) !=
                                             c("FAOST_CODE", "Year")],
                        collapse = " + "))
scatterplotMatrix(eval(parse(text = plotString)), data = imputed.dt,
                  smoother = FALSE, use = "pairwise.complete.obs")

## Merge with name table and set keys
imputed.dt = merge(imputed.dt, countryNamesTable, by = "FAOST_CODE")
imputed.dt = merge(imputed.dt,
    data.table(FAOregionProfile[, c("FAOST_CODE", "UNSD_MACRO_REG",
                                    "UNSD_SUB_REG")]), by = "FAOST_CODE")
setkeyv(imputed.dt, c("FAOST_CODE", "Year"))

## Post compute GDP per capita
imputed.dt[, BY.GDP.MKTP.PCAP :=
           log(exp(NY.GDP.MKTP.CD)/exp(POPULATION))]


## Merge country information to unimputed data for later comparison
final.dt = merge(final.dt, countryNamesTable, by = "FAOST_CODE")
final.dt = merge(final.dt,
    data.table(FAOregionProfile[, c("FAOST_CODE", "UNSD_MACRO_REG",
                                    "UNSD_SUB_REG")]), by = "FAOST_CODE")
final.dt[, BY.GDP.MKTP.PCAP := rep(NA, NROW(final.dt))]


## Back transform the data
## ---------------------------------------------------------------------

## logistic transformation for proportion
logit = function(x){
    1/(1 + exp(-x))
}


backTransformString = paste0("`:=`(",
    paste0(meta.dt[!is.na(TRANSFORM) & SOURCE == "raw", VARIABLE], " = ",
           ifelse(meta.dt[!is.na(TRANSFORM) & SOURCE == "raw",
                          TRANSFORM] == "logit", "logit(", "exp("),
           meta.dt[!is.na(TRANSFORM) & SOURCE == "raw", VARIABLE], ")",
           collapse = ", "), ")")
imputed.dt[, eval(parse(text = backTransformString))]

## ## These variables are not back transformed for visualization
## imputed.dt[, IS.RRS.DNST.K2 := logit(IS.RRS.DNST.K2)]
## imputed.dt[, IS.ROD.DNST.K2 := logit(IS.ROD.DNST.K2)]
## imputed.dt[, TI.IV.FEFTMT.USD.NO := exp(TI.IV.FEFTMT.USD.NO)]
## imputed.dt[, POPULATION := exp(POPULATION)]
## imputed.dt[, NY.GDP.MKTP.CD := exp(NY.GDP.MKTP.CD)]
## imputed.dt[, QV.NPV.FOOD.ID.NO := exp(QV.NPV.FOOD.ID.NO)]




## Q4:Does poverty redution always imply hunger reduction?
## ---------------------------------------------------------------------

xyplot(logit(POU) ~ logit(SI.POV.DDAY), data = final.dt,
       panel = function(x, y){
           panel.xyplot(x, y, type = c("g", "p", "r"))
           panel.abline(0, 1)
           }
       )


with(imputed.dt[Year == 2009, ],
     {
         plot(SI.POV.DDAY, POU, xlim = c(0, 1), ylim = c(0, 1),
              type = "n")
         text(SI.POV.DDAY, POU, labels = FAO_TABLE_NAME)
         abline(a = 0, b = 1, col = "red", lwd = 2)
         abline(a = 0, b = 1.5, col = "red", lty = 2)
         abline(a = 0, b = 1/1.5, col = "red", lty = 2)
     }
     )

with(imputed.dt[Year == 2009, ],
     {
         plot(SI.POV.NAHC, POU, xlim = c(0, 1), ylim = c(0, 1),
              type = "n")
         text(SI.POV.NAHC, POU, labels = FAO_TABLE_NAME)
         abline(a = 0, b = 1, col = "red", lwd = 2)
         abline(a = 0, b = 1.5, col = "red", lty = 2)
         abline(a = 0, b = 1/1.5, col = "red", lty = 2)
     }
     )

## Load and compare the nasometrix with world bank estimate
povNasometrix.df = read.csv("poverty_nasometrix.csv",
    stringsAsFactors = FALSE)
povNasometrix.df$level = with(povNasometrix.df,
    ifelse(right == "X", "right",
           ifelse(too.high == "X", "too.high",
                  ifelse(too.low == "X", "too.low", NA))))
povNasometrix.dt = data.table(fillCountryCode(data = povNasometrix.df,
    country = "country"))
povNasometrix.dt[country == "Cote d'Iviore", FAOST_CODE := 107]
povNasometrix.dt = povNasometrix.dt[!is.na(FAOST_CODE), ]
povNasometrix.dt[, Year := 2010]
povNasometrix.dt[, `:=`(right = NULL, too.high = NULL, too.low = NULL)]
setkeyv(povNasometrix.dt, c("FAOST_CODE", "Year"))
setnames(x = povNasometrix.dt, old = "should.be.about",
         new = "nasometrix")
povNasometrix.dt[, nasometrix := nasometrix/100]
povCompare.dt = merge(povNasometrix.dt,
    imputed.dt[, list(FAOST_CODE, Year, POU, SI.POV.DDAY, SI.POV.NAHC,
                      SI.POV.RUHC, SI.POV.GINI)], all.x = TRUE)

## Looks like the nasometrix are around the same with the world bank
## estimate, but adjust many countries which are below 5%
with(povCompare.dt,
     {
         plot(nasometrix, SI.POV.DDAY, xlim = c(0, 1), ylim = c(0, 1))
         abline(a = 0, b = 1)
     }
     )

## The nasometrix is almost stricly below the national poverty line
with(povCompare.dt,
     {
         plot(nasometrix, SI.POV.NAHC, xlim = c(0, 1), ylim = c(0, 1))
         abline(a = 0, b = 1)
     }
     )

nasometrixCountry = povCompare.dt[level == "right", FAOST_CODE]

## Keeping country where Josef and Piero consider that the poverty is
## valid.
with(imputed.dt[Year == 2009 & FAOST_CODE %in% nasometrixCountry, ],
     {
         plot(SI.POV.DDAY, POU, xlim = c(0, 1), ylim = c(0, 1),
              type = "n")
         text(SI.POV.DDAY, POU, labels = FAO_TABLE_NAME)
         abline(a = 0, b = 1, col = "red", lwd = 2)
         abline(a = 0, b = 1.5, col = "red", lty = 2)
         abline(a = 0, b = 1/1.5, col = "red", lty = 2)
     }
     )





## POU comparison between international and national poverty line
## pdf(file = "onetwentyfive_national_pou_comparison.pdf",
##     width = 20, height = 15)
par(mfrow = c(2, 1), mar = c(1, 4.1, 2.1, 2.1))
with(imputed.dt[Year == 2009, ],
     {
         plot(SI.POV.DDAY, POU, type = "n", xlim = c(0, 1),
              ylim = c(0, 1), axes = FALSE,
              ylab = "Prevalence of Undernourishment");
         box()
         text(SI.POV.DDAY, POU, labels = FAO_TABLE_NAME,
              cex = 0.8)
         text(1, 0.05, "Poverty headcount ratio at $1.25/day", adj = 1,
              cex = 1.5)
         abline(a = 0, b = 1, col = "red", lwd = 2)
         abline(a = 0, b = 1.5, col = "red", lty = 2)
         abline(a = 0, b = 1/1.5, col = "red", lty = 2)
     }
     )
legend("topleft",  legend = c("45 degeree line", "50% deviation"),
       lty = c(1, 2), lwd = c(1, 2), bty = "n", col = "red")
par(mar = c(3.1, 4.1, 0, 2.1))
with(imputed.dt[Year == 2009, ],
     {
         plot(SI.POV.NAHC, POU, type = "n", xlim = c(0, 1),
              ylim = c(0, 1),
              ylab = "Prevalence of Undernourishment");
         text(SI.POV.NAHC, POU, labels = FAO_TABLE_NAME,
              cex = 0.8)
         text(1, 0.05,
              "Poverty headcount ratio at national poverty line",
              adj = 1, cex = 1.5)
         abline(a = 0, b = 1, col = "red", lwd = 2)
         abline(a = 0, b = 1.5, col = "red", lty = 2)
         abline(a = 0, b = 1/1.5, col = "red", lty = 2)
     }
     )
## graphics.off()

## international versus national and rural poverty line
## pdf(file = "onetwentyfive_national_rural_comparison.pdf",
##     width = 20, height = 15)
par(mfrow = c(2, 1), mar = c(1, 4.1, 2.1, 2.1))
with(imputed.dt[Year == 2009, ],
     {
         plot(SI.POV.NAHC, SI.POV.DDAY, type = "n", xlim = c(0, 1),
              ylim = c(0, 1), axes = FALSE,
              ylab = "Poverty headcount ratio at $1.25/day");
         box()
         text(SI.POV.NAHC, SI.POV.DDAY, labels = FAO_TABLE_NAME,
              cex = 0.8)
         text(1, 0.05, "National poverty line", adj = 1, cex = 1.5)
         abline(a = 0, b = 1, col = "red", lwd = 2)
         abline(a = 0, b = 1.5, col = "red", lty = 2)
         abline(a = 0, b = 1/1.5, col = "red", lty = 2)
     }
     )
legend("topleft",  legend = c("45 degeree line", "50% deviation"),
       lty = c(1, 2), lwd = c(1, 2), bty = "n", col = "red")
par(mar = c(3.1, 4.1, 0, 2.1))
with(imputed.dt[Year == 2009, ],
     {
         plot(SI.POV.RUHC, SI.POV.DDAY, type = "n", xlim = c(0, 1),
              ylim = c(0, 1),
              ylab = "Poverty headcount ratio at $1.25/day");
         text(SI.POV.RUHC, SI.POV.DDAY, labels = FAO_TABLE_NAME,
              cex = 0.8)
         text(1, 0.05, "Rural poverty line", adj = 1, cex = 1.5)
         abline(a = 0, b = 1, col = "red", lwd = 2)
         abline(a = 0, b = 1.5, col = "red", lty = 2)
         abline(a = 0, b = 1/1.5, col = "red", lty = 2)
     }
     )
## graphics.off()


## examine countries in the lower left corner by log transform the data
with(imputed.dt[Year == 2009, ],
     {
         plot(log(SI.POV.NAHC), log(SI.POV.DDAY), type = "n",
              ylab = "Poverty headcount ratio at $1.25/day",
              xlab = "Poverty headcount ratio at national poverty line",
              axes = FALSE);
         axis(1, at = seq(-4, 0, by = 1),
              labels = round(exp(seq(-4, 0, by = 1)), 2))
         axis(2, at = seq(-15, 0, by = 5),
              labels = round(exp(seq(-15, 0, by = 5)), 5))
         box()
         text(log(SI.POV.NAHC), log(SI.POV.DDAY),
              labels = FAO_TABLE_NAME,
              cex = 0.8)
         text(1, 0.05, "National poverty line", adj = 1, cex = 1.5)
         abline(a = 0, b = 1, col = "red", lwd = 2)
         abline(a = 0, b = 1.5, col = "red", lty = 2)
         abline(a = 0, b = 1/1.5, col = "red", lty = 2)
     }
     )
legend("topleft",  legend = c("45 degeree line", "50% deviation"),
       lty = c(1, 2), lwd = c(1, 2), bty = "n", col = "red")

with(final.dt[Year %in% 2005:2009, ],
     {
         plot(log(logit(SI.POV.NAHC)), log(logit(SI.POV.DDAY)),
              type = "n",
              ylab = "Poverty headcount ratio at $1.25/day",
              xlab = "Poverty headcount ratio at national poverty line",
              axes = FALSE);
         axis(1, at = seq(-4, 0, by = 1),
              labels = round(exp(seq(-4, 0, by = 1)), 2))
         axis(2, at = seq(-15, 0, by = 5),
              labels = round(exp(seq(-15, 0, by = 5)), 5))
         box()
         text(log(logit(SI.POV.NAHC)), log(logit(SI.POV.DDAY)),
              labels = FAO_TABLE_NAME,
              cex = 0.8)
         text(1, 0.05, "National poverty line", adj = 1, cex = 1.5)
         abline(a = 0, b = 1, col = "red", lwd = 2)
         abline(a = 0, b = 1.5, col = "red", lty = 2)
         abline(a = 0, b = 1/1.5, col = "red", lty = 2)
     }
     )
legend("topleft",  legend = c("45 degeree line", "50% deviation"),
       lty = c(1, 2), lwd = c(1, 2), bty = "n", col = "red")



## This is a verification that the relationship is not a result of
## imputation
par(mfrow = c(2, 1), mar = c(1, 4.1, 2.1, 2.1))
with(final.dt[Year %in% 2005:2009, ],
     {
         plot(logit(SI.POV.NAHC), logit(SI.POV.DDAY), type = "n",
              xlim = c(0, 1), ylim = c(0, 1), axes = FALSE,
              ylab = "Poverty headcount ratio at $1.25/day");
         box()
         text(logit(SI.POV.NAHC), logit(SI.POV.DDAY),
              labels = FAO_TABLE_NAME, cex = 0.8)
         text(1, 0.05, "National poverty line", adj = 1, cex = 1.5)
         abline(a = 0, b = 1, col = "red", lwd = 2)
         abline(a = 0, b = 1.5, col = "red", lty = 2)
         abline(a = 0, b = 1/1.5, col = "red", lty = 2)
     }
     )
legend("topleft",  legend = c("45 degeree line", "50% deviation"),
       lty = c(1, 2), lwd = c(1, 2), bty = "n", col = "red")
par(mar = c(3.1, 4.1, 0, 2.1))
with(final.dt[Year %in% 2005:2009, ],
     {
         plot(logit(SI.POV.RUHC), logit(SI.POV.DDAY), type = "n",
              xlim = c(0, 1), ylim = c(0, 1),
              ylab = "Poverty headcount ratio at $1.25/day");
         text(logit(SI.POV.RUHC), logit(SI.POV.DDAY),
              labels = FAO_TABLE_NAME, cex = 0.8)
         text(1, 0.05, "Rural poverty line", adj = 1, cex = 1.5)
         abline(a = 0, b = 1, col = "red", lwd = 2)
         abline(a = 0, b = 1.5, col = "red", lty = 2)
         abline(a = 0, b = 1/1.5, col = "red", lty = 2)
     }
     )

## NOTE (Michael): Re-write this section so that the density
##                 comparison are specific to those year where
##                 undernourishment is more prevalent than poverty.




## Peculiar countries to be investigated
highPOU.dt = final.dt[logit(POU)/logit(SI.POV.DDAY) >= 2 &
    Year >= 2000 & logit(POU) >= 0.1 & logit(SI.POV.DDAY) >= 0.05,
    list(FAOST_CODE, FAO_TABLE_NAME, Year, logit(POU),
         logit(SI.POV.DDAY))]

unique(highPOU.dt[, list(FAOST_CODE, FAO_TABLE_NAME)])

## Examine country
examineVars = names(which(sapply(imputed.dt, typeof) == "double"))

x1 = rep(seq(0, 0.8, length = 6), 5)
y1 = rep(seq(0, 0.8, length = 5), each = 6)
x2 = rep(seq(0.2, 1, length = 6), 5)
y2 = rep(seq(0.2, 1, length = 5), each = 6)

n.vars = length(examineVars)



varNameTable = data.frame(examineVars = examineVars)
varNameTable.df = merge(varNameTable,
    meta.df[, c("STS_ID", "TITLE_STS")], by.x = "examineVars",
    by.y = "STS_ID", all.x = TRUE)
varNameTable.df[varNameTable.df$examineVars == "NY.GDP.MKTP.CD",
                "TITLE_STS"] = "GDP"
varNameTable.df[varNameTable.df$examineVars == "SI.POV.DDAY",
                "TITLE_STS"] = "Poverty headcount ratio at $1.25 per day"
varNameTable.df[varNameTable.df$examineVars == "SI.POV.GINI",
                "TITLE_STS"] = "Gini Index"
varNameTable.df[varNameTable.df$examineVars == "SI.POV.RUHC",
                "TITLE_STS"] =
    "Poverty headcount ratio at rural poverty line"
varNameTable.df[varNameTable.df$examineVars == "SI.POV.NAHC",
                "TITLE_STS"] =
    "Poverty headcount ratio at national poverty line"
varNameTable.df[varNameTable.df$examineVars == "SP.RUR.TOTL.ZS",
                "TITLE_STS"] = "Rural population percentage"
varNameTable.df[varNameTable.df$examineVars == "EN.POP.DNST",
                "TITLE_STS"] = "Population density"
## varNameTable.df[varNameTable.df$examineVars == "POP.ROD.DNST",
##                 "TITLE_STS"] = "Population to rode density"
varNameTable.df[varNameTable.df$examineVars == "BY.GDP.MKTP.PCAP",
                "TITLE_STS"] = "GDP per capita"

examineCountry = unique(imputed.dt$FAOST_CODE)
## examineCountry = unique(imputed.dt[IS.ROD.DNST.K2 <= -6 & IS.RRS.DNST.K2 <= -6 & ADESA <= 100 & SP.RUR.TOTL.ZS >= 0.7 & BY.GDP.MKTP.PCAP >= 6.5, FAOST_CODE])
examineCountry = unique(imputed.dt[Year == 2009 & SI.POV.DDAY < 0.1 & POU >= 0.1 & POU/SI.POV.DDAY >= 1.5, FAOST_CODE])
for(j in examineCountry){
    print(examineCountry)
    pdf(file = paste0("./exploratory_analysis/",
            FAOcountryProfile[which(FAOcountryProfile$FAOST_CODE == j),
                              "FAO_TABLE_NAME"], ".pdf"),
        width = 21 ,height = 21)
    axisPos.df = data.frame(x1 = x1, y1 = y1, x2 = x2, y2 = y2)
    for(i in seq_along(examineVars)){
        lab = varNameTable.df[varNameTable.df$examineVars ==
            examineVars[i], "TITLE_STS"]
        print(densityplot(eval(parse(text =
                                     paste0("~", examineVars[i]))),
                          data = imputed.dt[Year == 2010, ],
                          groups = ifelse(FAOST_CODE == j,
                              1, 0),
                          cex = c(1, 2),
                          main = list(label = lab, cex = 0.75),
                          xlab = "", ylab = ""),
              pos = axisPos.df[i, ], more = i != length(examineVars))
    }
    for(i in seq_along(examineVars)){
        lab = varNameTable.df[varNameTable.df$examineVars ==
            examineVars[i], "TITLE_STS"]
        print(xyplot(eval(parse(text =
                                paste0(examineVars[i], " ~ Year"))),
                     data = imputed.dt[FAOST_CODE == j, ],
                     type = c("g", "l"),
                     main = list(label = lab, cex = 0.75),
                     ylab = "", xlab = "",
                     ylim = range(c(imputed.dt[, examineVars[i],
                         with = FALSE]), na.rm = TRUE, finite = TRUE)),
              pos = axisPos.df[i, ], more = i != length(examineVars) + 1)
          print(xyplot(eval(parse(text =
                                paste0(examineVars[i], " ~ Year"))),
                       data = final.dt[FAOST_CODE == j, ],
                       type = c("g", "p"),
                       main = list(label = lab, cex = 0.75),
                       ylab = "", xlab = "",
                       ylim = range(c(imputed.dt[, examineVars[i],
                           with = FALSE]), na.rm = TRUE, finite = TRUE),
                       cex = 0.5),
                pos = axisPos.df[i, ], more = i != length(examineVars))
        if(examineVars[i] == "SI.POV.DDAY"){
          print(xyplot(logit(POU) + logit(SI.POV.DDAY) ~ Year,
                       data = final.dt[FAOST_CODE == j, ],
                       type = c("g", "p"),
                       main = list(label = lab, cex = 0.75),
                       ylab = "", xlab = "",
                       ylim = range(c(imputed.dt[, examineVars[i],
                           with = FALSE]), na.rm = TRUE, finite = TRUE),
                       cex = 0.5),
                pos = axisPos.df[i, ], more = TRUE)
      }
    }
    graphics.off()
}



## Old questions to be answered
## ---------------------------------------------------------------------



## Define categorical variables for interpretation
imputed.dt[, levelADESA := ifelse(ADESA < 90, "veryLow",
                            ifelse(ADESA >= 80 & ADESA < 100, "low",
                                   ifelse(ADESA >= 100 & ADESA < 115,
                                          "sufficient", "high")))]
imputed.dt[, levelADESA :=
           factor(levelADESA,
                  levels = c("veryLow", "low", "sufficient", "high"))]
imputed.dt[, levelSH.STA.STNT.ZS := ifelse(SH.STA.STNT.ZS < 0.05,
                                     "veryLow",
                            ifelse(SH.STA.STNT.ZS >= 0.05 &
                                   SH.STA.STNT.ZS < 0.1, "low",
                                   ifelse(SH.STA.STNT.ZS >= 0.1 &
                                          SH.STA.STNT.ZS < 0.3,
                                          "high", "veryHigh")))]
imputed.dt[, levelSH.STA.STNT.ZS :=
           factor(levelSH.STA.STNT.ZS,
                  levels = c("veryLow", "low", "high", "veryHigh"))]
imputed.dt[, levelSH.STA.STNT.ZS :=
           ifelse(SH.STA.STNT.ZS < 0.05, "veryLow",
                            ifelse(SH.STA.STNT.ZS >= 0.05 &
                                   SH.STA.STNT.ZS < 0.1, "low",
                                   ifelse(SH.STA.STNT.ZS >= 0.1 &
                                          SH.STA.STNT.ZS < 0.3,
                                          "high", "veryHigh")))]
imputed.dt[, levelSH.STA.STNT.ZS :=
           factor(levelSH.STA.STNT.ZS,
                  levels = c("veryLow", "low", "high", "veryHigh"))]
imputed.dt[, levelSH.STA.MALN.ZS :=
           ifelse(SH.STA.MALN.ZS < 0.02, "veryLow",
                            ifelse(SH.STA.MALN.ZS >= 0.02 &
                                   SH.STA.MALN.ZS < 0.05, "low",
                                   ifelse(SH.STA.MALN.ZS >= 0.05 &
                                          SH.STA.MALN.ZS < 0.2,
                                          "high", "veryHigh")))]
imputed.dt[, levelSH.STA.MALN.ZS :=
           factor(levelSH.STA.MALN.ZS,
                  levels = c("veryLow", "low", "high", "veryHigh"))]
imputed.dt[, levelSH.STA.WAST.ZS :=
           ifelse(SH.STA.WAST.ZS < 0.025, "veryLow",
                            ifelse(SH.STA.WAST.ZS >= 0.02 &
                                   SH.STA.WAST.ZS < 0.05, "low",
                                   ifelse(SH.STA.WAST.ZS >= 0.05 &
                                          SH.STA.WAST.ZS < 0.15,
                                          "high", "veryHigh")))]
imputed.dt[, levelSH.STA.WAST.ZS :=
           factor(levelSH.STA.WAST.ZS,
                  levels = c("veryLow", "low", "high", "veryHigh"))]



## Q1:Does improved access to food also mean better utilization?
## ---------------------------------------------------------------------
##
## A1: Short answer is yes, but it depends on the level of food
##     availability and current level of stunting.


ggplot(imputed.dt, aes(x = ADESA, y = SH.STA.STNT.ZS)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "lm", formula = y ~ log(x)) +
    facet_wrap(~Year)


xyplot(POPULATION ~Year|FAO_TABLE_NAME, data = imputed.dt,
       type = c("g", "l"))


xyplot(SH.STA.STNT.ZS ~ ADESA, data = imputed.dt,
       type = c("g", "p", "smooth"))

xyplot(SH.STA.STNT.ZS + ADESA/100 ~Year|FAO_TABLE_NAME, data = imputed.dt,
       type = c("g", "l"), auto.key = TRUE)

xyplot(POU ~ SH.STA.STNT.ZS, data = imputed.dt)

imputed.dt[, diffADESA := c(NA, diff(ADESA)), by = "FAOST_CODE"]
imputed.dt[, diffSH.STA.STNT.ZS := c(NA, diff(SH.STA.STNT.ZS)),
           by = "FAOST_CODE"]
imputed.dt[, diffSH.STA.WAST.ZS := c(NA, diff(SH.STA.WAST.ZS)),
           by = "FAOST_CODE"]
imputed.dt[, diffSH.STA.MALN.ZS := c(NA, diff(SH.STA.MALN.ZS)),
           by = "FAOST_CODE"]
## imputed.dt[, diffSH.STA.AMALN.ZS := c(NA, diff(SH.STA.AMALN.ZS)),
##            by = "FAOST_CODE"]



lagpad <- function(x, k) {
    c(rep(NA, k), x)[1 : length(x)]
}
imputed.dt[, diffLag1ADESA := lagpad(diffADESA, 1),
           by = "FAOST_CODE"]
imputed.dt[, diffLag3ADESA := lagpad(diffADESA, 3),
           by = "FAOST_CODE"]
imputed.dt[, diffLag5ADESA := lagpad(diffADESA, 5),
           by = "FAOST_CODE"]



foo = function(formula){
    tmp = try(coef(lm(formula))[2])
    if(inherits(tmp, "try-error"))
        tmp = NA
    as.numeric(tmp)
}

imputed.dt[, regressionSlope := foo(diffSH.STA.STNT.ZS ~ diffADESA),
           by = c("levelADESA", "levelSH.STA.STNT.ZS")]


ggplot(data = imputed.dt[!is.na(levelADESA) & !is.na(levelSH.STA.STNT.ZS), ],
       aes(x = diffADESA, y = diffSH.STA.STNT.ZS)) +
    geom_point(alpha = 0.1) +
    stat_smooth(method="lm", se = FALSE, fullrange = TRUE) +
    facet_grid(levelADESA ~ levelSH.STA.STNT.ZS) +
    geom_text(aes(x = 0, y = 0.75,
                  label = paste0("Slope = ", ifelse(is.na(regressionSlope), 0,
                      round(regressionSlope, 3)))))

ggplot(data = imputed.dt[!is.na(levelADESA) & !is.na(levelSH.STA.STNT.ZS) &
           abs(diffADESA) < 20 & abs(diffSH.STA.STNT.ZS) < 20, ],
       aes(x = diffADESA, y = diffSH.STA.STNT.ZS)) +
    geom_point(alpha = 0.1) +
    stat_smooth(method="lm", se = FALSE, fullrange = TRUE) +
    facet_grid(levelADESA ~ levelSH.STA.STNT.ZS)


ggplot(data = imputed.dt[!is.na(levelADESA) & !is.na(levelSH.STA.WAST.ZS) &
           abs(diffADESA) < 20 & abs(diffSH.STA.WAST.ZS) < 20, ],
       aes(x = diffADESA, y = diffSH.STA.WAST.ZS)) +
    geom_point(alpha = 0.1) +
    stat_smooth(method="lm", se = FALSE, fullrange = TRUE) +
    facet_grid(levelADESA ~ levelSH.STA.WAST.ZS)


ggplot(data = imputed.dt[!is.na(levelADESA) & !is.na(levelSH.STA.MALN.ZS), ],
       aes(x = diffADESA, y = diffSH.STA.MALN.ZS)) +
    geom_point(alpha = 0.1) +
    stat_smooth(method="lm", se = FALSE, fullrange = TRUE) +
    facet_grid(levelADESA ~ levelSH.STA.MALN.ZS)

ggplot(data = imputed.dt[!is.na(levelADESA) & !is.na(levelSH.STA.MALN.ZS) &
           abs(diffADESA) < 20 & abs(diffSH.STA.MALN.ZS) < 20, ],
       aes(x = diffADESA, y = diffSH.STA.MALN.ZS)) +
    geom_point(alpha = 0.1) +
    stat_smooth(method="lm", se = FALSE, fullrange = TRUE) +
    facet_grid(levelADESA ~ levelSH.STA.MALN.ZS)




## Q2:Does high food availability impy lower undernourishment?
## ---------------------------------------------------------------------
##
## A2: This is by construction, we don't actually know the real rate
##     of undernourishment

xyplot(POU ~ ADESA|FAO_TABLE_NAME, data = imputed.dt, type = c("g", "p", "r"))

xyplot(POU ~ ADESA|FAO_TABLE_NAME, data = final.dt, type = c("g", "p", "r"))


## What is the relationship between accessability and undernourishment?
## ---------------------------------------------------------------------
##
## Infrastructure paves the way for reducing undernourishment
xyplot(POU ~ IS.ROD.DNST.K2, data = imputed.dt)
xyplot(POU ~ IS.ROD.DNST.K2|FAO_TABLE_NAME, data = imputed.dt)
xyplot(logit(POU) ~ IS.ROD.DNST.K2|FAO_TABLE_NAME, data = final.dt)



ggplot(data = imputed.dt[!is.na(levelADESA), ],
       aes(x = IS.ROD.DNST.K2, y = POU)) +
    geom_point() + facet_grid(Year ~levelADESA)


ggplot(data = imputed.dt[!is.na(levelADESA) & IS.ROD.DNST.K2 <= 10, ],
       aes(x = IS.ROD.DNST.K2, y = POU)) +
    geom_point() + facet_wrap(~Year)



xyplot(log(POU) + log(IS.ROD.DNST.K2)~Year|FAO_TABLE_NAME,
       data = imputed.dt[IS.ROD.DNST.K2 <= 10, ],
       type = c("g", "l"), auto.key = TRUE)

xyplot(logit(POU) ~ IS.ROD.DNST.K2|levelADESA,
       data = imputed.dt,
       xlim = c(0, 10))

xyplot(logit(POU) ~ IS.ROD.DNST.K2|Year, data = final.dt,
       xlim = c(0, 10))


xyplot(POU ~ IS.ROD.DNST.K2|levelADESA, data = imputed.dt,
       xlim = c(0, 10))

xyplot(POU ~ IS.ROD.DNST.K2|FAO_TABLE_NAME, data = imputed.dt,
       xlim = c(0, 10))



xyplot(IS.RRS.DNST.K2 ~ IS.ROD.DNST.K2|FAO_TABLE_NAME, data = imputed.dt,
       xlim = c(0, 10))


xyplot(POU ~ IS.RRS.DNST.K2, data = imputed.dt)
xyplot(POU ~ IS.RRS.DNST.K2|FAO_TABLE_NAME, data = imputed.dt)

xyplot(logit(POU) + IS.RRS.DNST.K2 ~ Year|FAO_TABLE_NAME, data = final.dt)



xyplot(IS.ROD.DNST.K2 + IS.RRS.DNST.K2 ~ Year|FAO_TABLE_NAME, data = final.dt,
       type = c("g", "l"), ylim = c(0, 10))







xyplot(SI.POV.DDAY + SI.POV.RUHC ~ Year|FAO_TABLE_NAME,
       data = imputed.dt, type = c("g", "l"))



xyplot(logit(SI.POV.DDAY) + logit(SI.POV.RUHC) + logit(SI.POV.NAHC) ~
       Year|FAO_TABLE_NAME.x,
       data = final.dt, type = c("g", "l"), auto.key = TRUE)

xyplot(logit(SI.POV.DDAY) ~ logit(SI.POV.RUHC),
       data = final.dt, type = c("g", "p"), auto.key = TRUE)


xyplot(POU ~ SI.POV.DDAY|FAO_TABLE_NAME, data = imputed.dt,
       panel = function(x, y){
           panel.xyplot(x, y, type = c("g", "p", "r"), cex = 0.5)
           panel.abline(0, 1)
       },
       type = c("g", "p", "r"), auto.key = TRUE, xlim = c(0, 1),
       ylim = c(0, 1), cex = 0.5)


xyplot(POU + SI.POV.DDAY ~ Year|FAO_TABLE_NAME, data = imputed.dt,
       type = c("g", "l"), auto.key = TRUE)


xyplot(POU - SI.POV.DDAY ~ Year|FAO_TABLE_NAME, data = imputed.dt,
       type = c("g", "l"))


xyplot(logit(POU) + logit(SI.POV.DDAY) ~ Year|FAO_TABLE_NAME,
       data = final.dt,
       type = c("g", "l"), auto.key = TRUE)



























q1.dt = final.dt[ , list(FAOST_CODE, FAO_TABLE_NAME, UNSD_SUB_REG,
    UNSD_MACRO_REG, Year, ADESA, POU, SH.STA.WAST.ZS, SH.STA.MALN.ZS,
    SH.STA.STNT.ZS)]
popU14 = getWDItoSYB("SP.POP.0014.TO.ZS")
q1Final.dt = data.table(mergeSYB(data.frame(q1.dt),
    popU14$entity[, c("ISO2_WB_CODE", "Year", "SP.POP.0014.TO.ZS")]))
q1Final.dt = q1Final.dt[Year >= 1985, ]

q1Final.dt = test3

xyplot(SP.POP.0014.TO.ZS ~ Year|FAO_TABLE_NAME, data = q1Final.dt, type = "l")

q1Final.dt[, ageClass := ifelse(SP.POP.0014.TO.ZS < 23, "low",
                          ifelse(SP.POP.0014.TO.ZS > 42, "high",
                                 "medium"))]


with(q1Final.dt, plot(POU * 100, SH.STA.STNT.ZS,
                      col = rgb(ifelse(is.na(SP.POP.0014.TO.ZS), 0,
                          SP.POP.0014.TO.ZS/100), 0, 0, alpha = 0.75),
                      pch = 19))

## The old plot
xyplot(POU ~ SH.STA.STNT.ZS|FAO_TABLE_NAME, data = q1Final.dt,
       type = c("p", "r"))

xyplot(POU * 100 ~ SH.STA.STNT.ZS|ageClass, data = q1Final.dt,
       type = c("p", "r"), xlim = c(0, 100), ylim = c(0, 100))


xyplot(POU + SH.STA.STNT.ZS ~ Year|FAO_TABLE_NAME, data = q1Final.dt,
       type = c("l"), auto.key = TRUE)

xyplot(SH.STA.WAST.ZS + SH.STA.MALN.ZS + SH.STA.STNT.ZS ~
       Year|FAO_TABLE_NAME, data = q1Final.dt,
       type = c("l"), auto.key = TRUE)

xyplot(diffADESA ~ diffSH.STA.STNT.ZS|FAO_TABLE_NAME, data = q1Final.dt)


## Imputed stunt
## q1Final.dt[, STNTPOURATIO := SH.STA.STNT.ZS/(POU * 100)]
## q1Final.dt[STNTPOURATIO >= 10, STNTPOURATIO := as.numeric(NA)]
## stuntImputationModel = lmer(STNTPOURATIO ~
##     (-1 + Year|FAO_TABLE_NAME),
##     data = q1Final.dt)

stuntImputationModel = lmer(SH.STA.STNT.ZS ~
    (-1 + POU + log(ADESA)|FAO_TABLE_NAME),
    data = q1Final.dt)

## q1Final.dt[is.na(STNTPOURATIO),
##            STNTPOURATIO := predict(stuntImputationModel, newdata = .SD,
##                             allow.new.levels = TRUE)]
## q1Final.dt[!is.na(SH.STA.STNT.IMP),
##            SH.STA.STNT.IMP := STNTPOURATIO * (POU * 100)]

q1Final.dt[!is.na(POU) & !is.na(ADESA),
           SH.STA.STNT.IMP := predict(stuntImputationModel, newdata = .SD,
                               allow.new.levels = TRUE)]
xyplot(POU ~ SH.STA.STNT.IMP|as.factor(Year), data = q1Final.dt)
xyplot(ADESA ~ SH.STA.STNT.IMP|as.factor(Year), data = q1Final.dt)

pctdiff = function(x){
    c(NA, diff(x))/x[length(x)]
}

q1Final.dt[, diffPOU := c(NA, diff(POU)), by = "FAO_TABLE_NAME"]
q1Final.dt[, diffADESA := c(NA, diff(ADESA)), by = "FAO_TABLE_NAME"]
q1Final.dt[, diffSH.STA.STNT.ZS := c(NA, diff(SH.STA.STNT.ZS)),
           by = "FAO_TABLE_NAME"]
q1Final.dt[, diffLag1SH.STA.STNT.ZS :=
           lagpad(diffSH.STA.STNT.ZS, k = 1), by = "FAO_TABLE_NAME"]
q1Final.dt[, diffLag3SH.STA.STNT.ZS :=
           lagpad(diffSH.STA.STNT.ZS, k = 3), by = "FAO_TABLE_NAME"]
q1Final.dt[, diffLag5SH.STA.STNT.ZS :=
           lagpad(diffSH.STA.STNT.ZS, k = 5), by = "FAO_TABLE_NAME"]

q1Final.dt[, diffSH.STA.STNT.ZS := c(NA, diff(SH.STA.STNT.ZS)),
           by = "FAO_TABLE_NAME"]
q1Final.dt[, sufficientADESA := ifelse(ADESA >= 100, "Yes", "No")]
q1Final.dt[, levelADESA := ifelse(ADESA < 80, "veryLow",
                            ifelse(ADESA >= 80 & ADESA < 100, "low",
                                   ifelse(ADESA >= 100 & ADESA < 115,
                                          "sufficient", "high")))]
q1Final.dt[, levelADESA :=
           factor(levelADESA,
                  levels = c("veryLow", "low", "sufficient", "high"))]

q1Final.dt[, levelSH.STA.STNT.ZS := ifelse(SH.STA.STNT.ZS < 0.05, "veryLow",
                            ifelse(SH.STA.STNT.ZS >= 0.05 &
                                   SH.STA.STNT.ZS < 0.15, "low",
                                   ifelse(SH.STA.STNT.ZS >= 0.15 &
                                          SH.STA.STNT.ZS < 0.25,
                                          "high", "veryHigh")))]
q1Final.dt[, levelSH.STA.STNT.ZS :=
           factor(levelSH.STA.STNT.ZS,
                  levels = c("veryLow", "low", "high", "veryHigh"))]


xyplot(diffPOU ~ diffSH.STA.STNT.IMP|sufficientADESA, data = q1Final.dt)
xyplot(diffADESA ~ diffSH.STA.STNT.IMP|sufficientADESA, data = q1Final.dt)
xyplot(diffADESA ~ diffLag1SH.STA.STNT.IMP|sufficientADESA, data = q1Final.dt)
xyplot(diffADESA ~ diffLag3SH.STA.STNT.IMP|sufficientADESA, data = q1Final.dt)
xyplot(diffADESA ~ diffLag5SH.STA.STNT.IMP|sufficientADESA, data = q1Final.dt)
xyplot(diffADESA ~ diffSH.STA.STNT.ZS|sufficientADESA, data = q1Final.dt)

xyplot(ADESA + SH.STA.STNT.IMP ~ Year|FAO_TABLE_NAME, data = q1Final.dt,
       type = c("l", "g"), auto.key = TRUE)
xyplot(diffADESA ~ diffSH.STA.STNT.IMP|FAO_TABLE_NAME, data = q1Final.dt)


xyplot(diffADESA ~ diffSH.STA.STNT.IMP|sufficientADESA, data = q1Final.dt)
xyplot(diffADESA ~ diffSH.STA.STNT.IMP|levelADESA, data = q1Final.dt)
xyplot(diffADESA ~ diffLag1SH.STA.STNT.IMP|levelADESA, data = q1Final.dt)
xyplot(diffADESA ~ diffLag3SH.STA.STNT.IMP|levelADESA, data = q1Final.dt)
xyplot(diffADESA ~ diffLag5SH.STA.STNT.IMP|levelADESA, data = q1Final.dt)

xyplot(diffADESA ~ diffPOU|levelADESA, data = q1Final.dt, type = c("p", "r"))
xyplot(diffADESA ~ diffSH.STA.STNT.IMP|levelADESA, data = q1Final.dt,
       type = c("p", "r"))



print(xyplot(diffSH.STA.STNT.ZS ~ diffADESA|levelADESA,
       data = q1Final.dt[SH.STA.STNT.ZS >= 0.05, ],
       type = c("g", "p", "r")))

print(xyplot(SH.STA.STNT.ZS ~ ADESA,
       data = q1Final.dt[SH.STA.STNT.ZS >= 0.05, ],
       type = c("g", "p", "r")))


print(xyplot(diffSH.STA.STNT.ZS ~ diffADESA|levelSH.STA.STNT.ZS,
       data = q1Final.dt[SH.STA.STNT.ZS >= 0.05, ],
       type = c("g", "p", "r")))


print(xyplot(diffSH.STA.STNT.ZS ~ diffADESA|levelSH.STA.STNT.ZS,
       data = q1Final.dt[SH.STA.STNT.ZS >= 0.05, ],
       type = c("g", "p", "r")))



ggplot(data = q1Final.dt[SH.STA.STNT.ZS >= 0.05 & !is.na(levelADESA) &
           diffADESA >= -40,],
       aes(x = diffADESA, y = diffSH.STA.STNT.ZS)) +
    geom_point(alpha = 0.2) +
    stat_smooth(method="lm", se = FALSE, fullrange = TRUE) +
    facet_grid(levelADESA ~ levelSH.STA.STNT.ZS) +
    scale_y_continuous(limits=c(-50, 50)) +
    scale_x_continuous(limits=c(-50, 50))



ggplot(data = q1Final.dt[SH.STA.STNT.ZS >= 5 & !is.na(levelADESA) &
           diffADESA >= -40,],
       aes(x = diffADESA, y = diffSH.STA.STNT.ZS, col = levelSH.STA.STNT.ZS)) +
    geom_point() + stat_smooth(method="lm", se = FALSE, fullrange = TRUE) +
    facet_grid(~ levelADESA) +
    scale_y_continuous(limits=c(-50, 50)) +
    scale_x_continuous(limits=c(-50, 50))


print(xyplot(ADESA ~ SH.STA.STNT.IMP|levelADESA,
       data = q1Final.dt[SH.STA.STNT.IMP >= 0, ],
       type = c("g", "l", "r")))


print(xyplot(ADESA ~ SH.STA.STNT.IMP,
       data = q1Final.dt[SH.STA.STNT.IMP >= 0, ],
       type = c("g", "l"), group = FAO_TABLE_NAME))


with(q1Final.dt, plot(diffADESA, diffSH.STA.STNT.IMP))


xyplot(POU * 100 + SH.STA.STNT.IMP ~ Year|FAO_TABLE_NAME, data = q1Final.dt,
       type = c("l"), auto.key = TRUE)

xyplot(POU * 100 ~ SH.STA.STNT.IMP|FAO_TABLE_NAME, data = q1Final.dt,
       type = c("p", "r"), xlim = c(0, 100), ylim = c(0, 100))



xyplot(POU * 100 + SH.STA.STNT.IMP ~ Year|FAO_TABLE_NAME,
       data = q1Final.dt[FAO_TABLE_NAME != "Kuwait" & POU < 0.1, ],
       type = c("l"), auto.key = TRUE)


xyplot(ADESA/3 + SH.STA.STNT.IMP ~ Year|FAO_TABLE_NAME,
       data = q1Final.dt[FAO_TABLE_NAME != "Kuwait" & POU < 0.2, ],
       type = c("l"), auto.key = TRUE)



xyplot(STNTPOURATIO ~ Year|FAO_TABLE_NAME, data = q1Final.dt,
       type = c("g", "l"))
