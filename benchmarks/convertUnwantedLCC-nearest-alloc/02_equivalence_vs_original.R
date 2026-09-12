## Validate the EDITED in-package convertUnwantedLCC vs the INSTALLED original.
suppressPackageStartupMessages({ library(terra); library(data.table); library(LandR) })
paddedFloatToChar <- reproducible::paddedFloatToChar
CLONE <- Sys.getenv("LANDR_SRC", ".")

## ---- parse check + eval the edited function ----
ex <- parse(file.path(CLONE, "R/cohorts.R"))            # errors here if syntax broken
cuMod <- NULL
for (e in ex) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                  identical(e[[2]], as.name("convertUnwantedLCC"))) { cuMod <- eval(e); break }
stopifnot(is.function(cuMod)); cat("edited convertUnwantedLCC parsed + evaluated OK\n")

applyOut <- function(r,o){ v<-as.vector(r[]); o<-o[!is.na(pixelIndex)]; v[o$pixelIndex]<-suppressWarnings(as.integer(o$ecoregionGroup)); setValues(rast(r),v) }
match_uw <- function(a,b,r,unw){ uw<-which(as.vector(r[])==unw); 100*mean(as.vector(a[])[uw]==as.vector(b[])[uw],na.rm=TRUE) }

## ===== TEST 1: unconstrained (plain LCC codes) =====
n<-90; r<-rast(nrows=n,ncols=n,xmin=0,xmax=n,ymin=0,ymax=n); xy<-xyFromCell(r,1:ncell(r))
values(r)<-ifelse(xy[,1]<n/2,210L,220L); d<-sqrt((xy[,1]-n/2)^2+(xy[,2]-n/2)^2); r[d<25]<-240L
set.seed(2); r[sample(ncell(r),120)]<-240L
aDT<-data.table(pixelIndex=seq_len(ncell(r)), initialEcoregionCode=as.integer(values(r)[,1]))
set.seed(10); oOld1<-suppressMessages(LandR::convertUnwantedLCC(240L,r,copy(aDT),doAssertion=FALSE))
set.seed(11); oOld2<-suppressMessages(LandR::convertUnwantedLCC(240L,r,copy(aDT),doAssertion=FALSE))
oMod <-suppressMessages(cuMod(240L,r,copy(aDT),doAssertion=FALSE))
oMod2<-suppressMessages(cuMod(240L,r,copy(aDT),doAssertion=FALSE))
cat(sprintf("TEST1 unconstrained: old-self=%.1f%%  mod-vs-old=%.1f%%  mod resolves all=%s  deterministic=%s\n",
  match_uw(applyOut(r,oOld1),applyOut(r,oOld2),r,240), match_uw(applyOut(r,oMod),applyOut(r,oOld1),r,240),
  sum(values(applyOut(r,oMod))[,1]==240,na.rm=TRUE)==0, identical(oMod,oMod2)))

## ===== TEST 2: constrained (per-ecoregion availability, preDash codes) =====
eco<-ifelse(xy[,1]<n/2,"1","2"); lcc<-ifelse(xy[,1]<n/2,210L,220L); dd<-sqrt((xy[,1]-n/2)^2+(xy[,2]-n/2)^2)
lcc[dd<25]<-240L; set.seed(4); lcc[sample(ncell(r),120)]<-240L
rC<-setValues(rast(r),lcc); iec<-paste0(eco,"_",formatC(lcc,width=3,flag="0"))
aDT2<-data.table(pixelIndex=seq_len(ncell(rC)), initialEcoregionCode=iec)
set.seed(20); oOldC1<-suppressMessages(LandR::convertUnwantedLCC(240L,rC,copy(aDT2),doAssertion=FALSE))
set.seed(21); oOldC2<-suppressMessages(LandR::convertUnwantedLCC(240L,rC,copy(aDT2),doAssertion=FALSE))
oModC<-suppressMessages(cuMod(240L,rC,copy(aDT2),doAssertion=FALSE))
cmp <-merge(oOldC1[!is.na(pixelIndex)],oModC,by="pixelIndex",suffixes=c(".old",".mod"))
cmpOO<-merge(oOldC1[!is.na(pixelIndex)],oOldC2[!is.na(pixelIndex)],by="pixelIndex",suffixes=c(".1",".2"))
cat(sprintf("TEST2 constrained  : old-self=%.1f%%  mod-vs-old=%.1f%%  mod never assigns unavailable ERC=%s\n",
  100*mean(cmpOO$ecoregionGroup.1==cmpOO$ecoregionGroup.2), 100*mean(cmp$ecoregionGroup.old==cmp$ecoregionGroup.mod),
  all(oModC[!is.na(ecoregionGroup)]$ecoregionGroup %in% aDT2$initialEcoregionCode)))
cat("assigned-ERC distribution (old vs mod):\n")
print(merge(oOldC1[!is.na(ecoregionGroup),.N,ecoregionGroup], oModC[!is.na(ecoregionGroup),.N,ecoregionGroup],
            by="ecoregionGroup",all=TRUE,suffixes=c(".old",".mod")))
