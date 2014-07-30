##  loosely derived from library(goseq); example(goseq)

gseseq <- function(genes, isDE, weights, sets)
    ## genes: character(n)
    ## isDE: logical(n)
    ## weights: numeric(n); all(0 <= weights)
    ## sets: list of gene names
{
    ## validate (a little) inputs
    inp <- dta.frame(isDE=isDE, weight=weights, row.names=genes)
    stopifnot(all(unlist(sets, use.names=FALSE) %in% rownames(inp)))

    ## overall summaries
    nTot <- nrow(inp)
    nTotDE <- sum(inp$isDE)
    alpha <- sum(inp$weight)

    p <- sapply(sets, function(elt, nTot, nTotDE) {
        ## gene set calculations
        nInSet <- length(elt)
        nInSetDE <- sum(inp[elt, "isDE"])
        avewt <- mean(inp[elt, "weight"])
        wt <- (nTot - nInSet) * avewt / (alpha - nInSet * avewt)
        ## Wallenius tests
        d <- dWNCHypergeo(nInSetDE, nInSet, nTot - nInSet, nTotDE, wt)
        p0 <- pWNCHypergeo(nInSetDE, nInSet, nTot - nInSet, nTotDE, wt,
                           lower.tail = FALSE)
        p1 <- pWNCHypergeo(nInSetDE, nInSet, nTot - nInSet, nTotDE, wt)
        ## return
        c(over=d + p0, under=p1)        
    }, nTot, nTotDE)

    as.data.frame(t(p))
}
