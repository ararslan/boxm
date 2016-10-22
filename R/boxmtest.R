boxmtest <- function(X, groups) {
    if (!inherits(X, "matrix") && !inherits(X, "data.frame")) {
        stop("X must be a numeric matrix or data frame")
    }
    if (!is.factor(groups)) {
        stop("groups must be a factor")
    }
    if (length(groups) != nrow(X)) {
        stop("groups must have length equal to the number of rows in X")
    }
    if (nlevels(groups) < 2) {
        stop("groups has fewer than 2 levels")
    }

    N <- NROW(X)
    p <- NCOL(X)

    n <- table(groups)
    nl <- as.list(setNames(n, levels(groups)))

    g <- nlevels(groups)
    Ng <- N - g

    Si <- lapply(levels(groups), function(i) cov(X[groups == i, ]))
    Sp <- Reduce("+", lapply(1:g, function(i) (nl[[i]] - 1) * Si[[i]])) / Ng

    u <- (sum(1 / (n - 1)) - 1 / Ng) * ((2 * p^2 + 3 * p - 1) / (6 * (p + 1) * (g - 1)))
    M <- Ng * log(det(Sp)) - sum((n - 1) * sapply(Si, function(x) log(det(x))))

    df <- as.integer((p * (p + 1) * (g - 1)) %/% 2)

    if (p <= 5 && g <= 5) {
        stat <- (1 - u) * M
        return(list(statistic = stat, df = df, pvalue = pchisq(stat, df, lower.tail = TRUE)))
    } else {
        a <- (sum(1 / (n - 1)^2) - 1 / Ng^2) * (p^2 + p - 2) / (6 * (g - 1))
        df2 <- as.integer((df + 2) / abs(a - u^2))

        if (a > u^2) {  # Pearson's type VI curve
            b <- df / (1 - u - df / df2)
            stat <- M / b
        } else {        # Pearson's type I curve
            b <- df / (1 - u - 2 / df2)
            stat <- df2 * M / (df * (b - M))
        }
        # Hopefully the type III case (a == u^2) will be taken care of by the chisq case
        # above, since df2 will be infinite so the statistic is ~ chisq(df)

        return(list(stat = stat, df = c(df, df2), pvalue = pf(stat, df, df2, lower.tail = TRUE)))
    }
}
