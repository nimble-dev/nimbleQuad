## TODO: these might be moved to `nimble`.

makeNodeString <- function(nodes, model) {
    if (!length(nodes))
        return("")
    elements <- model$expandNodeNames(nodes, returnScalarComponents = TRUE)
    vars <- sapply(strsplit(elements, "[", fixed = TRUE), `[[`, 1)
    nodesCount <- table(vars)
    start <- names(nodesCount)
    end <- c(rep(", ", length(start) - 1), "")
    end[nodesCount > 1] <- paste0(" (", nodesCount[nodesCount > 1], " elements)",
        end[nodesCount > 1])
    return(paste0(c(rbind(start, end)), collapse = ""))
}

splitLatents <- function(model, paramNodes, latentNodes, calcNodes, calcNodesOther,
    control = list()) {
    stochNodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
    discreteStochNodes <- model$isDiscrete(stochNodes)
    if (any(discreteStochNodes))
        stop("splitLatents: found discrete non-data stochastic nodes in processing nodes for quadrature-based posterior approximation: ",
            paste0(stochNodes[discreteStochNodes], collapse = ", "), ". Discrete non-data stochastic nodes cannot be handled by the posterior approximation algorithm.")
    split <- extractControlElement(control, "split", TRUE)
    check <- extractControlElement(control, "check", TRUE)
    margNodes <- setupMargNodes(model = model, paramNodes = paramNodes, randomEffectsNodes = latentNodes,
        calcNodes = calcNodes, calcNodesOther = calcNodesOther, split = split, check = check)
    if (missing(paramNodes) && missing(latentNodes)) {
        if (!missing(calcNodes) || !missing(calcNodesOther))
            messageIfVerbose("   [Note] Ignoring provide `calcNodes` and `calcNodesOther` because `paramNodes` and `latentNodes` not provided and are being determinted automatically.")
        paramNodes <- margNodes$paramNodes
        latentNodes <- margNodes$randomEffectsNodes
        deps <- model$getDependencies(latentNodes, includeData = FALSE, self = FALSE)
        ## By default, we treat "siblings" of latent nodes as latents.
        ## This attempts to have fixed effects in latents,
        ## along with random effects.
        newLatents <- model$getParents(deps)
        paramNodes <- setdiff(paramNodes, newLatents)
        latentNodes <- unique(c(latentNodes, newLatents))
        margNodes <- setupMargNodes(model = model, paramNodes = paramNodes,
            randomEffectsNodes = latentNodes, split = split, check = check)
    }

    return(margNodes)
}

nodesToIndices <- function(node, approx) {
    Rapprox <- ifelse(is(approx, "NestedApprox"), approx, approx$Robject)
    if (!node %in% Rapprox$paramNodesComponents)
        stop("`", node, "` not found in model")
    index <- paramNodesIndices[which(node == paramNodesComponents)]
    if (index == 0)
        stop("`", node, "` is part of a multivariate node for which the parameter transformation is not one-to one. Use `sampleLatentNodes()` to sample from the approximate posterior for the parameters.")
    return(index)
}
