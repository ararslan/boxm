boxmtest <- function(X, groups) {
    if (!inherits(X, "matrix") && !inherits(X, "data.frame")) {
        stop("X must be a numeric matrix or data frame")
    }
    if (!is.factor(groups)) {
        stop("groups must be a factor.")
    }
    if (length(groups) != nrow(X)) {
        stop("groups must have length equal to the number of rows in X")
    }

    N <- NROW(X)
    p <- NCOL(X)

    n <- table(groups)
    nl <- as.list(setNames(n, levels(groups)))

    g <- nlevels(groups)

    Si <- lapply(levels(groups), function(i) cov(X[groups == i, ]))
    Sp <- Reduce("+", lapply(1:g, function(i) (nl[[i]] - 1) * Si[[i]])) / N

    u <- sum(1 / (n - 1) - 1 / N) * ((2 * p^2 + 3 * p - 1) / (6 * (p + 1) * (g - 1)))
    M <- N * log(det(Sp)) - sum((n - 1) * sapply(Si, function(x) log(det(x))))

    stat <- (1 - u) * M
    df <- p * (p + 1) * (g - 1) / 2

    list(statistic = stat, df = df, pvalue = pchisq(stat, df))
}
