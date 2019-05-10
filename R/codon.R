#' converts codon nucleotide sequences to amino acids
#' 
#' @param x first sequence
#' @param y default = NULL. If provided, will print aa1>aa2.
#' @return converted amino acid annotation
#' @examples
#' codon("ttt")
#' # [1] "Phe"
#' codon("ttt", "ttg")
#' # [1] "Phe>Leu"
#' @export

codon <- function(x, y = NULL)
{
    if (toupper(x) == "TTT" | toupper(x) == "TTC")
    {
        res <- "Phe"
    }
    if (toupper(x) == "TTA" | toupper(x) == "TTG" | toupper(x) == "CTT" | 
        toupper(x) == "CTC" | toupper(x) == "CTA" | toupper(x) == "CTG")
        {
        res <- "Leu"
    }
    if (toupper(x) == "ATT" | toupper(x) == "ATC" | toupper(x) == "ATA")
    {
        res <- "Ile"
    }
    if (toupper(x) == "ATG")
    {
        res <- "Met"
    }
    if (toupper(x) == "GTT" | toupper(x) == "GTC" | toupper(x) == "GTA" | 
        toupper(x) == "GTG")
        {
        res <- "Val"
    }
    if (toupper(x) == "TCT" | toupper(x) == "TCC" | toupper(x) == "TCA" | 
        toupper(x) == "TCG")
        {
        res <- "Ser"
    }
    if (toupper(x) == "CCT" | toupper(x) == "CCC" | toupper(x) == "CCA" | 
        toupper(x) == "CCG")
        {
        res <- "Pro"
    }
    if (toupper(x) == "ACT" | toupper(x) == "ACC" | toupper(x) == "ACA" | 
        toupper(x) == "ACG")
        {
        res <- "Thr"
    }
    if (toupper(x) == "GCT" | toupper(x) == "GCC" | toupper(x) == "GCA" | 
        toupper(x) == "GCG")
        {
        res <- "Ala"
    }
    if (toupper(x) == "TAT" | toupper(x) == "TAC")
    {
        res <- "Tyr"
    }
    if (toupper(x) == "TAA" | toupper(x) == "TAG" | toupper(x) == "TGA")
    {
        res <- "Stop"
    }
    if (toupper(x) == "CAT" | toupper(x) == "CAC")
    {
        res <- "His"
    }
    if (toupper(x) == "CAA" | toupper(x) == "CAG")
    {
        res <- "Gln"
    }
    if (toupper(x) == "AAT" | toupper(x) == "AAC")
    {
        res <- "Asn"
    }
    if (toupper(x) == "AAA" | toupper(x) == "AAG")
    {
        res <- "Lys"
    }
    if (toupper(x) == "GAT" | toupper(x) == "GAC")
    {
        res <- "Asp"
    }
    if (toupper(x) == "GAA" | toupper(x) == "GAG")
    {
        res <- "Glu"
    }
    if (toupper(x) == "TGT" | toupper(x) == "TGC")
    {
        res <- "Cys"
    }
    if (toupper(x) == "TGG")
    {
        res <- "Trp"
    }
    if (toupper(x) == "CGT" | toupper(x) == "CGC" | toupper(x) == "CGA" | 
        toupper(x) == "CGG" | toupper(x) == "AGA" | toupper(x) == "AGG")
        {
        res <- "Arg"
    }
    if (toupper(x) == "AGT" | toupper(x) == "AGC")
    {
        res <- "Ser"
    }
    if (toupper(x) == "GGT" | toupper(x) == "GGC" | toupper(x) == "GGA" | 
        toupper(x) == "GGG")
        {
        res <- "Gly"
    }
    
    if (length(y) != 0)
    {
        if (toupper(y) == "TTT" | toupper(y) == "TTC")
        {
            res2 <- "Phe"
        }
        if (toupper(y) == "TTA" | toupper(y) == "TTG" | toupper(y) == "CTT" | 
            toupper(y) == "CTC" | toupper(y) == "CTA" | toupper(y) == "CTG")
            {
            res2 <- "Leu"
        }
        if (toupper(y) == "ATT" | toupper(y) == "ATC" | toupper(y) == "ATA")
        {
            res2 <- "Ile"
        }
        if (toupper(y) == "ATG")
        {
            res2 <- "Met"
        }
        if (toupper(y) == "GTT" | toupper(y) == "GTC" | toupper(y) == "GTA" | 
            toupper(y) == "GTG")
            {
            res2 <- "Val"
        }
        if (toupper(y) == "TCT" | toupper(y) == "TCC" | toupper(y) == "TCA" | 
            toupper(y) == "TCG")
            {
            res2 <- "Ser"
        }
        if (toupper(y) == "CCT" | toupper(y) == "CCC" | toupper(y) == "CCA" | 
            toupper(y) == "CCG")
            {
            res2 <- "Pro"
        }
        if (toupper(y) == "ACT" | toupper(y) == "ACC" | toupper(y) == "ACA" | 
            toupper(y) == "ACG")
            {
            res2 <- "Thr"
        }
        if (toupper(y) == "GCT" | toupper(y) == "GCC" | toupper(y) == "GCA" | 
            toupper(y) == "GCG")
            {
            res2 <- "Ala"
        }
        if (toupper(y) == "TAT" | toupper(y) == "TAC")
        {
            res2 <- "Tyr"
        }
        if (toupper(y) == "TAA" | toupper(y) == "TAG" | toupper(y) == "TGA")
        {
            res2 <- "Stop"
        }
        if (toupper(y) == "CAT" | toupper(y) == "CAC")
        {
            res2 <- "His"
        }
        if (toupper(y) == "CAA" | toupper(y) == "CAG")
        {
            res2 <- "Gln"
        }
        if (toupper(y) == "AAT" | toupper(y) == "AAC")
        {
            res2 <- "Asn"
        }
        if (toupper(y) == "AAA" | toupper(y) == "AAG")
        {
            res2 <- "Lys"
        }
        if (toupper(y) == "GAT" | toupper(y) == "GAC")
        {
            res2 <- "Asp"
        }
        if (toupper(y) == "GAA" | toupper(y) == "GAG")
        {
            res2 <- "Glu"
        }
        if (toupper(y) == "TGT" | toupper(y) == "TGC")
        {
            res2 <- "Cys"
        }
        if (toupper(y) == "TGG")
        {
            res2 <- "Trp"
        }
        if (toupper(y) == "CGT" | toupper(y) == "CGC" | toupper(y) == "CGA" | 
            toupper(y) == "CGG" | toupper(y) == "AGA" | toupper(y) == "AGG")
            {
            res2 <- "Arg"
        }
        if (toupper(y) == "AGT" | toupper(y) == "AGC")
        {
            res2 <- "Ser"
        }
        if (toupper(y) == "GGT" | toupper(y) == "GGC" | toupper(y) == "GGA" | 
            toupper(y) == "GGG")
            {
            res2 <- "Gly"
        }
        out <- paste0(res, ">", res2)
    } else
    {
        out <- paste0(res)
    }
    kelvinny::pbcopy(out)
    return(out)
}