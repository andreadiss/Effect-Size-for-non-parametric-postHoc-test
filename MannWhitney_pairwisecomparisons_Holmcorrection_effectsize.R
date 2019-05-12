# -------------------------------------------------------------------------------------------------------------------------------------------------------#
# In this script I changed the original pairwise.wilcox.test function to compute two indices of the effect size for the Wilcox test (and Mann-Whitney U)
# for repeated comparisons. The two indices are based on the Probability Superiority and the U standardization procedures. General reference to the two procedures
# can be find in the book: Effect sizes for research by Krissom and Kim.

# The first two tables represent the U statistic and the p-value - Holm corrected - per comparisons.
# The third table represent the the denominator of the Probability of Superiority index (namely na*nb).
# The fourth table is the Probability of superiority index.
# The fifth table is the r index computed after standardizing the U statistics.
#
# The arguments of the pairwise2.wilcox.test are: x = DV, g = group variable, n = sample size
# -------------------------------------------------------------------------------------------------------------------------------------------------------#


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
