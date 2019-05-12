# Questo script calcola nelle prime due tavole il p.value (Holm's correction) e la statistica W
# la terza tavola è il prodotto n1*n2, la quarta rappresenta la probabilità di superiorità (PS) e la quinta tavola l'effect size per la statistica U

# inserire la VD, la variabile gruppo e N del campione

# DCSF:per analisi
# library(PMCMRplus)
# dscfAllPairsTest(x, g)

pairwise2.table<-function (funct, level.names)
{
    ix <- setNames(seq_along(level.names), level.names)
    pp <- outer(ix[-1L], ix[-length(ix)], function(ivec, jvec) sapply(seq_along(ivec), 
        function(k) {
            i <- ivec[k]
            j <- jvec[k]
            if (i > j)  funct(i, j, n) else NA
        }))
    pp
}



pairwise2.wilcox.test<-function (x, g, n, p.adjust.method = "holm", paired = FALSE) 
{
    p.adjust.method <- match.arg(p.adjust.method)
    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
    g <- factor(g)	
    METHOD <- if (paired)  "Wilcoxon signed rank test" else "Wilcoxon rank sum test"
    compare.levels <- function(i, j, n) {
        xi <- x[as.integer(g) == i]
        xj <- x[as.integer(g) == j]
        round(wilcox.test(xi, xj, paired = paired)$p.value,4)
    }  
    compare2.levels <- function(i, j, n) {
        xi <- x[as.integer(g) == i]
        xj <- x[as.integer(g) == j]
        wilcox.test(xi, xj, paired = paired)$statistic
    }
   compare3.levels <- function(i, j, n) {
        xi <- x[as.integer(g) == i]
        xj <- x[as.integer(g) == j]
        length(xi)*length(xj)
    }

   compare4.levels <- function(i, j, n) {
        xi <- x[as.integer(g) == i]
        xj <- x[as.integer(g) == j]
        wilcox.test(xi, xj, paired = paired)$statistic/(length(xi)*length(xj))
    }

   compare5.levels <- function(i, j, n) {
        xi <- x[as.integer(g) == i]
        xj <- x[as.integer(g) == j]
        num = wilcox.test(xi, xj, paired = paired)$statistic-((length(xi)*length(xj))/2)
	  z = abs(num)/sqrt((length(xi)*length(xj)*(length(xi)+length(xj)+1))/12)
	  z/sqrt(n)
    }


    PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
    Uw<-pairwise2.table(funct=compare2.levels, levels(g))
    DEN<-pairwise2.table(funct=compare3.levels, levels(g))
    PS <- pairwise2.table(funct=compare4.levels, levels(g))
    EffectSize_method2 <- pairwise2.table(funct=compare5.levels, levels(g))
	
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL, W=Uw,DEN=DEN, PS = PS, EffectSize_method2=EffectSize_method2,
        p.adjust.method = p.adjust.method)
        ans
}
