## This file is for binomial distributions

#' #' EM algorithm for esitmating binomial mixture model
#' #' @param k A vector of integers. number of success
#' #' @param n A vector of integers. number of trials
#' #' @param n_components A number. number of components
#' #' @param tol numeric(1), tolerance value for convergence between iterations
#' #' @return p1 A float between 0 and 1. mean success fraction of the first binomial
#' #' @return p2 A float between 0 and 1. mean success fraction of the second binomial
#' #' @return psi A float between 0 and 1. fraction of the first component
#' #' @examples
#' #' n1 <- array(sample(1:30, 50, replace = TRUE))
#' #' n2 <- array(sample(1:30, 200, replace = TRUE))
#' #' k1 <- apply(n1, 1, rbinom, n = 1, p = 0.5)
#' #' k2 <- apply(n2, 1, rbinom, n = 1, p = 0.01)
#' #' RV <- mixBinom(c(k1, k2), c(n1, n2))
#' mixBinom <- function(k, n, n_components = 2) {
#'     S <- length(n)
#'     ## Random initialzation on parameters
#'     p <- runif(n_components, 0.0, 1.0)
#'     psi <- runif(n_components, 0.0, 1.0)
#'     psi <- psi / sum(psi)
#'
#'     p_new <- runif(n_components, 0.0, 1.0)
#'     psi_new <- runif(n_components, 0.0, 1.0)
#'     psi_new <- psi_new / sum(psi_new)
#'
#'     prob_mat <- matrix(runif(n_components * S, 0, 1), nrow = S, byrow = TRUE)
#'     prob_mat <- prob_mat / rowSums(prob_mat)
#'
#'     ## Iterations
#'     while (!(all(abs(p - p_new) < tol))) {
#'         ## E-step:
#'         for (j in seq_len(n_components)) {
#'                 prob_mat[,j] <- psi[j] * dbinom(k, size = n, prob = p_new[j],
#'                                                 log = FALSE)
#'         }
#'         prob_mat <- prob_mat / rowSums(prob_mat)
#'
#'         ## M-step
#'         for (j in seq_len(n_components)) {
#'             p[j] <- p_new[j]
#'             p_new[j] <- sum(prob_mat[,j] * k) / sum(prob_mat[,j] * n)
#'             psi[j] <- psi_new[j]
#'             psi_new[j] <- sum(prob_mat[,j]) / S
#'         }
#'     }
#'
#'     ## return values
#'     return_list <- list("p" = p_new, "psi" = psi_new, "prob" = prob_mat)
#'     return_list
#'  }
#'
#'
#'  #' Predicted probability from learned binomial mixture model
#'  #' @examples
#'  #' RV <- mixBinom(c(k1,k2), c(n1,n2))
#'  #' prob <- pred_mixBinom(3, 10, RV$p, RV$psi)
#'  pred_mixBinom <- function(k, n, p, psi){
#'    prob_test <- runif(length(p), 0, 1)
#'    for (j in 1:length(p)){
#'      prob_test[j] <- dbinom(k, size=n, prob=p[j], log = FALSE) * psi[j]
#'    }
#'    prob_test <- prob_test / sum(prob_test)
#'    prob_test
#'  }


get_clone_score <- function(x, var_meta, clone = "clone1") {
    locs <- dplyr::filter(var_meta, clone == clone)[["locName"]]
    score <- rowMeans(as.matrix(x[, locNames(x) %in% locs]) > 0, na.rm = TRUE)
    score[is.nan(score)] <- 0
    score[is.na(score)] <- 0
    score
}

compute_relatedness_vec <- function(geno_vec, geno_mat, vaf = NULL) {
    ## geno_mat should be samples x variants
    if ( !is.null(vaf) ) {
        ## standardise genotypes
        message("Standardise genotypes\n")
        stand_geno_vec <- ( (geno_vec - 2 * vaf)
            / sqrt(2 * vaf * (1 - vaf)) )
        stand_geno_mat<- t( (t(geno_mat) - 2 * vaf)
                                 / sqrt(2 * vaf * (1 - vaf)))
        message("Compute genotypic correlations\n")
        corr_mat <- t(t(stand_geno_mat) * stand_geno_vec)
        #bad_snp <- apply(corr_mat, 2, function(x) any(is.na(x)))
        #corr_mat <- corr_mat[, !bad_snp, drop = FALSE]
        ## sum correlations across donors
        ave_allelic_corr <- rowMeans(corr_mat, na.rm = TRUE)
    } else {
        ave_allelic_corr <- rep(NA, length(geno_vec))
    }
    ave_allelic_corr
}

simil_realrel <- function(geno_mat, vaf) {
    ## geno_mat should be samples x variants
    simil_out <- apply(geno_mat, 1, compute_relatedness_vec, geno_mat, vaf)
    simil_out
}



dist_shared_alleles_vec <- function(geno_vec, geno_mat) {
    ## compute distance between a sample and a set of samples
    ## based on alternative allele count
    raw_mat <-  abs(t(t(geno_mat) - geno_vec)) / 2
    ave_allelic_dist <- rowMeans(raw_mat, na.rm = TRUE)
}

dist_shared_alleles <- function(geno_mat) {
    dist_out <- apply(geno_mat, 1, dist_shared_alleles_vec, geno_mat)
    as.dist(dist_out)
}


plot_sil_width <- function(sill) {
    df <- data_frame(cluster = as.factor(sill[, 1]),
                     neighbor = as.factor(sill[, 2]),
                     sil_width = sill[, 3])
    df <- dplyr::arrange(df, cluster, sil_width) %>%
        dplyr::mutate(y = 1:nrow(.))
    ave_wid <- mean(df[["sil_width"]])
    df_ave <- group_by(df, cluster) %>%
        summarise(ave_sil_wid = mean(sil_width),
                  y_pos = mean(y)) %>%
        dplyr::mutate(
            lab = paste0("Ave. ", format(ave_sil_wid, digits = 2)))
    ggplot(df, aes(y = y, x = sil_width, colour = cluster)) +
        geom_segment(aes(xend = 0, yend = y), size = 0.5) +
        geom_text(aes(x = -0.7, y = y_pos, label = lab), data = df_ave,
                  show.legend = FALSE) +
        geom_vline(xintercept = 0, linetype = 2) +
        scale_color_tableau(palette = "tableau20") +
        xlab("Silhouette width") +
        theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
              axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
        ggtitle(paste0("Ave. silh. width: ", format(ave_wid, digits = 3)))
}

ggmds <- function(dist_mat, df, colour = "cluster") {
    mds_out <- cmdscale(dist_mat, k = 6)
    df <- bind_cols(
        df, data_frame(MDS_1 = mds_out[, 1], MDS_2 = mds_out[, 2],
                       MDS_3 = mds_out[, 3], MDS_4 = mds_out[, 4],
                       MDS_5 = mds_out[, 5], MDS_6 = mds_out[, 6]))
    p <- ggplot(df,
           aes_string(x = "MDS_1", y = "MDS_2", colour = colour)) +
        geom_point(alpha = 0.5) + theme(legend.position = "bottom") +
        guides(size = FALSE)
    print(p)
    df
}


gensimil <- function(dat, lambda = mean(dat > 0, na.rm = TRUE), delta = 1,
                     nan_replace = 0.001) {
    bindat <- (dat > 0)
    m11 <- (replace(bindat, is.na(bindat), 0) %*%
            t(replace(bindat, is.na(bindat), 0)))
    m00 <- (replace(!bindat, is.na(bindat), 0) %*%
            t(replace(!bindat, is.na(bindat), 0)))
    m01 <- (replace(!bindat, is.na(bindat), 0) %*%
            t(replace(bindat, is.na(bindat), 0)))
    m10 <- (replace(bindat, is.na(bindat), 0) %*%
            t(replace(!bindat, is.na(bindat), 0)))
    sim <- ( (m11 + lambda * m00) /
             (m11 + lambda * m00 + delta * (m10 + m01)) )
    sim[is.nan(sim)] <- nan_replace
    diag(sim) <- 1
    sim
}

get_ap_clusters <- function(ap) {
    out <- as.factor(ap@idx)
    for (i in seq_along(levels(out)))
         levels(out) <- gsub(levels(out)[i], i, levels(out))
    out
}


canopy_post <- function(sampchain, projectname, K, numchain, burnin, thin,
                        optK, C = NULL, post.config.cutoff = NULL) {
    if (is.null(C)) {
        C = diag(nrow(sampchain[[1]][[1]][[1]]$cna))
        colnames(C) = rownames(C) = rownames(sampchain[[1]][[1]][[1]]$cna)
    }
    if (is.null(post.config.cutoff)) {
        post.config.cutoff = 0.05
    }
    if (post.config.cutoff > 1 | post.config.cutoff <= 0) {
        stop("post.config.cutoff has to be between 0 and 1!")
    }
    sampchaink = sampchain[[which(K == optK)]]
    numchain = length(sampchaink)
    samptreenew = sampchaink[[1]][(burnin + 1):length(sampchaink[[1]])]
    numpostburn = length(samptreenew)
    temp <- thin * c(1:(numpostburn/thin))
    samptreethin = samptreenew[temp]
    length(samptreethin)
    for (numi in 2:numchain) {
        samptreenew = sampchaink[[numi]][(burnin + 1):length(sampchaink[[numi]])]
        numpostburn = length(samptreenew)
        temp <- thin * c(1:(numpostburn/thin))
        samptreethin = c(samptreethin, samptreenew[temp])
    }
    samptreethin.lik = rep(NA, length(samptreethin))
    for (treei in 1:length(samptreethin)) {
        samptreethin.lik[treei] = samptreethin[[treei]]$likelihood
    }
    samptreethin = samptreethin[which((rank(-1 * samptreethin.lik,
                                            ties.method = "first")) <= 5 * (length(samptreethin)/numchain))]
    samptreethin.lik = rep(NA, length(samptreethin))
    for (treei in 1:length(samptreethin)) {
        samptreethin.lik[treei] = samptreethin[[treei]]$likelihood
    }
    if (!is.null(sampchain[[1]][[1]][[1]]$cna)) {
        for (i in 1:length(samptreethin)) {
            samptreethin[[i]] = sortcna(samptreethin[[i]], C)
        }
    }
    for (i in 1:length(samptreethin)) {
        samptreethin[[i]]$clonalmut = getclonalcomposition(samptreethin[[i]])
    }
    config = rep(NA, length(samptreethin))
    config[1] = 1
    categ = 1
    for (i in 2:length(samptreethin)) {
        for (categi in 1:categ) {
            list.a = samptreethin[[i]]$clonalmut
            list.b = samptreethin[[which(config == categi)[1]]]$clonalmut
            if ((sum(is.element(list.a, list.b)) == optK) &
                (sum(is.element(list.b, list.a)) == optK)) {
                config[i] = categi
            }
        }
        if (is.na(config[i])) {
            config[i] = categ + 1
            categ = categ + 1
        }
    }
    z.temp = (samptreethin.lik - mean(samptreethin.lik))/sd(samptreethin.lik)
    samptreethin = samptreethin[z.temp <= 1.5 & z.temp >= -1.5]
    samptreethin.lik = samptreethin.lik[z.temp <= 1.5 & z.temp >=
                                        -1.5]
    config = config[z.temp <= 1.5 & z.temp >= -1.5]
    config.summary = matrix(nrow = length(unique(config)), ncol = 3)
    colnames(config.summary) = c("Configuration", "Post_prob",
                                 "Mean_post_lik")
    config.summary[, 1] = unique(config)
    for (i in 1:nrow(config.summary)) {
        configi = config.summary[i, 1]
        configi.temp = which(config == configi)
        config.summary[i, 2] = round(length(configi.temp)/length(config),
                                     3)
        config.summary[i, 3] = round(max(samptreethin.lik[which(config ==
                                                                configi)]), 2)
    }
    minor.config = which(config.summary[, 2] < post.config.cutoff)
    if (length(minor.config) > 0) {
        config.sel = rep(TRUE, length(config))
        for (i in minor.config) {
            config.sel[which(config == config.summary[i, 1])] = FALSE
        }
        samptreethin = samptreethin[config.sel]
        samptreethin.lik = samptreethin.lik[config.sel]
        config = config[config.sel]
        config.summary = config.summary[-minor.config, , drop = FALSE]
        for (i in 1:nrow(config.summary)) {
            config[which(config == config.summary[i, 1])] = i
        }
        config.summary[, 1] = 1:nrow(config.summary)
        config.summary[, 2] = round(config.summary[, 2]/sum(config.summary[,
                                                                           2]), 3)
    }
    for (treei in 1:length(samptreethin)) {
        output.tree = samptreethin[[treei]]
        output.tree.Z = output.tree$Z[, 2:ncol(output.tree$Z), drop = FALSE]
        output.tree.P = apply(
            output.tree$P[2:nrow(output.tree$P),, drop = FALSE], 2,
            function(x) { x/sum(x) })
        output.tree$CCF = output.tree.Z %*% output.tree.P
        samptreethin[[treei]] = output.tree
    }
    return(list(samptreethin, samptreethin.lik, config, config.summary))
}


canopy_output_to_df <- function(tree) {
    p <- data_frame(clone = rownames(tree$P), prevalence = tree$P[, 1])
    outdf <- as_data_frame(tree$Z) %>%
        dplyr::mutate(var_id = rownames(tree$Z), vaf = tree$VAF[, 1],
                      ccf = tree$CCF[, 1])
    outdf[["clone"]] <- colnames(tree$Z)[1]
    for (i in 2:ncol(tree$Z))
        outdf[["clone"]][as.logical(tree$Z[, i])] <- colnames(tree$Z)[i]
    outdf
}

# plot_tree <- function(tree, title = NULL, size = 10) {
#     ptree <- ggtree(tree, ladderize = FALSE, layout = "slanted") +
#         geom_tiplab(size = 10, color = "firebrick")
#     df_prev <- as.data.frame(tree$P)
#     df_prev[["clone"]] <- rownames(tree$P)
#     df_prev[["branch.y"]] <- ptree$data$branch.y[1:sum(!is.na(ptree$data$label))]
#     df_prev[["branch.y.f"]] <- as.factor(as.numeric(as.factor(df_prev[["branch.y"]])))
#     df_prev <- gather(df_prev, key = "sample", value = "prevalence", -clone,
#                       -branch.y, -branch.y.f)
#     df_prev <- dplyr::mutate(df_prev, prevalence = signif(prevalence, digits = 3))
#     ptile <- ggplot(df_prev, aes(x = sample, y = branch.y.f)) +
#         geom_tile(aes(fill = prevalence)) +
#         geom_label(aes(label = prevalence), colour = "firebrick", size = 5) +
#         scale_fill_viridis(option = "B") +
#         scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
#         theme_minimal() + theme(axis.title = element_blank(),
#                                 axis.text.y = element_blank())
#     pt <- plot_grid(ptree, ptile, nrow = 1, rel_widths = c(1, 0.5))
#     if (!is.null(title)) {
#         title <- ggdraw() + draw_label(title, fontface='bold')
#         pt <- plot_grid(title, pt, ncol = 1, rel_heights = c(0.1, 1))
#     }
#     pt
# }


################################################################################

initialsna.bugfix <- function(tree, sna.name) {
    sna.no = length(sna.name)
    sna.edge = sample(2:nrow(tree$edge), size = sna.no, replace = TRUE)
    sna.mat = cbind(1:sna.no, tree$edge[sna.edge,, drop=FALSE])
    colnames(sna.mat) = c("sna", "sna.st.node", "sna.ed.node")
    rownames(sna.mat) = sna.name
    return(sna.mat)
}

canopy.sample.cluster.nocna.bugfix = function(R, X, sna_cluster, K, numchain,
                                       max.simrun, min.simrun, writeskip,
                                       projectname, cell.line = NULL,
                                       plot.likelihood = NULL) {
    if(length(sna_cluster)!=nrow(R)){
        stop('Length of sna_cluster should be the same as row numbers of R and X!')
    }
    if (!is.matrix(R)) {
        stop("R should be a matrix!")
    }
    if (!is.matrix(X)) {
        stop("X should be a matrix!")
    }
    if (min(K) < 2) {
        stop("Smallest number of subclones should be >= 2!\n")
    }
    if (is.null(cell.line)) {
        cell.line = FALSE
    }
    if (is.null(plot.likelihood)) {
        plot.likelihood = TRUE
    }
    if ( plot.likelihood){
        pdf(file = paste(projectname, "_likelihood.pdf", sep = ""), width = 10, height = 5)
    }

    sampname = colnames(R)
    sna.name = rownames(R)
    sampchain = vector("list", length(K))
    ki = 1
    for (k in K) {
        cat("Sample in tree space with", k, "subclones\n")
        sampchaink = vector("list", numchain)
        sampchaink.lik=vector('list',numchain)
        sampchaink.accept.rate=vector('list',numchain)
        for (numi in 1:numchain) {  # numi: number of chain
            cat("\tRunning chain", numi, "out of", numchain, "...\n")
            ###################################### Tree initialization #####
            text = paste(paste(paste(paste("(", 1:(k - 1), ",", sep = ""),
                                     collapse = ""), k, sep = ""), paste(rep(")", (k - 1)),
                                                                         collapse = ""), ";", sep = "")
            runif.temp=runif(1)
            if(k == 5 & runif.temp<0.5){
                text = c('(1,((2,3),(4,5)));')
            }else if(k == 6 & runif.temp < 1/3){
                text = c('(1,((2,3),(4,(5,6))));')
            }else if(k == 6 & runif.temp > 2/3){
                text = c('(1,(2,((3,4),(5,6))));')
            }else if(k == 7 & runif.temp > 1/4 & runif.temp <= 2/4){
                text=c('(1,((2,3),(4,(5,(6,7)))));')
            }else if(k == 7 & runif.temp > 2/4 & runif.temp <= 3/4){
                text = c('(1,((2,3),((4,5),(6,7))));')
            }else if(k == 7 & runif.temp > 3/4){
                text = c('(1,((2,(3,4)),(5,(6,7))));')
            }
            tree <- read.tree(text = text)
            tree$sna.cluster=initialsna.bugfix(tree,
                paste('cluster',sort(unique(sna_cluster)),sep=''))
            sna.mat = cbind(sna=1:nrow(R),
                (tree$sna.cluster)[paste0("cluster", sna_cluster),2:3])
            colnames(sna.mat) = c("sna", "sna.st.node", "sna.ed.node")
            rownames(sna.mat) = sna.name
            tree$sna=sna.mat
            #tree$sna = initialsna(tree, sna.name)
            # if(k>=5){tree$relation=getrelation(tree)}
            tree$Z = getZ(tree, sna.name)
            tree$P = initialP(tree, sampname, cell.line)
            tree$VAF = tree$Z%*%tree$P/2
            tree$likelihood = getlikelihood.sna(tree, R, X)
            ###################################### Sample in tree space #####
            sampi = 1
            writei = 1
            samptree = vector("list", max.simrun)
            samptree.lik=rep(NA, max.simrun)
            samptree.accept=rep(NA, max.simrun)
            samptree.accept.rate=rep(NA, max.simrun)

            while(sampi <= min.simrun){
                ######### sample sna mutation cluster positions
                tree.new=tree
                tree.new$sna.cluster=sampsna.cluster(tree)
                sna.mat = cbind(sna=1:nrow(R),(tree.new$sna.cluster)[paste0("cluster", sna_cluster),2:3])
                colnames(sna.mat) = c("sna", "sna.st.node", "sna.ed.node")
                rownames(sna.mat) = sna.name
                tree.new$sna=sna.mat
                tree.new$Z = getZ(tree.new, sna.name)
                tree.new$VAF=tree.new$Z%*%tree.new$P/2
                tree.new$likelihood = getlikelihood.sna(tree.new, R, X)
                tree.temp=addsamptree(tree,tree.new)
                tree=tree.temp[[1]]
                samptree.accept[sampi]=tree.temp[[2]]
                if (sampi%%writeskip == 0) {
                    samptree[[writei]] = tree
                    writei = writei + 1
                }
                samptree.lik[sampi]=tree$likelihood
                samptree.accept.rate[sampi]=mean(samptree.accept[max(1,sampi-999):sampi])
                sampi = sampi + 1
                ######## sample P (clonal proportions)
                tree.new = tree
                tree.new$P = sampP(tree.new, cell.line)
                tree.new$VAF = tree.new$Z%*%tree.new$P/2
                tree.new$likelihood = getlikelihood.sna(tree.new, R, X)
                tree.temp=addsamptree(tree,tree.new)
                tree=tree.temp[[1]]
                samptree.accept[sampi]=tree.temp[[2]]
                if (sampi%%writeskip == 0) {
                    samptree[[writei]] = tree
                    writei = writei + 1
                }
                samptree.lik[sampi]=tree$likelihood
                samptree.accept.rate[sampi]=mean(samptree.accept[max(1,sampi-999):sampi])
                sampi = sampi + 1
            }
            while(sampi <= max.simrun){
                ######### sample sna positions
                tree.new = tree
                tree.new$sna = sampsna(tree)
                tree.new$Z = getZ(tree.new, sna.name)
                tree.new$VAF=tree.new$Z%*%tree.new$P/2
                tree.new$likelihood = getlikelihood.sna(tree.new, R, X)
                tree.temp=addsamptree(tree,tree.new)
                tree=tree.temp[[1]]
                samptree.accept[sampi]=tree.temp[[2]]
                if (sampi%%writeskip == 0) {
                    samptree[[writei]] = tree
                    writei = writei + 1
                }
                samptree.lik[sampi]=tree$likelihood
                samptree.accept.rate[sampi]=mean(samptree.accept[max(1,sampi-999):sampi])
                if ((sampi >= 2*min.simrun) & (samptree.lik[sampi] <= mean(samptree.lik[max((sampi-1000),1):max((sampi-1),1)])) &
                    (samptree.accept.rate[sampi] <= mean(samptree.accept.rate[max((sampi-1000),1):max((sampi-1),1)])) &
                    (samptree.accept.rate[sampi] <= 0.1)) break
                sampi = sampi + 1
                ######## sample P (clonal proportions)
                tree.new = tree
                tree.new$P = sampP(tree.new, cell.line)
                tree.new$VAF = tree.new$Z%*%tree.new$P/2
                tree.new$likelihood = getlikelihood.sna(tree.new, R, X)
                tree.temp=addsamptree(tree,tree.new)
                tree=tree.temp[[1]]
                samptree.accept[sampi]=tree.temp[[2]]
                if (sampi%%writeskip == 0) {
                    samptree[[writei]] = tree
                    writei = writei + 1
                }
                samptree.lik[sampi]=tree$likelihood
                samptree.accept.rate[sampi]=mean(samptree.accept[max(1,sampi-999):sampi])
                if ((sampi >= 2*min.simrun) & (samptree.lik[sampi] <= mean(samptree.lik[max((sampi-1000),1):max((sampi-1),1)])) &
                    (samptree.accept.rate[sampi] <= mean(samptree.accept.rate[max((sampi-1000),1):max((sampi-1),1)])) &
                    (samptree.accept.rate[sampi] <= 0.1)) break
                sampi = sampi + 1
            }
            sampchaink[[numi]] = samptree[1:(writei - 1)]
            sampchaink.lik[[numi]]=samptree.lik
            sampchaink.accept.rate[[numi]]=samptree.accept.rate
        }
        ###################################### plotting and saving #####
        if (plot.likelihood) {
            par(mfrow=c(1,2))
            xmax=ymin=ymax=rep(NA,numchain)
            for(i in 1:numchain){
                xmax[i]=max(which((!is.na(sampchaink.lik[[i]]))))
                ymin[i]=sampchaink.lik[[i]][1]
                ymax[i]=sampchaink.lik[[i]][xmax[i]]
            }
            plot(sampchaink.lik[[1]],xlim=c(1,max(xmax)),ylim=c(min(ymin),max(ymax)),
                 xlab='Iteration',ylab='Log-likelihood',type='l',
                 main=paste('Post. likelihood:',k,'branches'))
            for(numi in 2:numchain){
                points(sampchaink.lik[[numi]],xlim=c(1,max(xmax)),ylim=c(min(ymin),max(ymax)),col=numi,type='l')
            }
            plot(sampchaink.accept.rate[[1]],ylim=c(0,1),xlim=c(1,max(xmax)),
                 xlab='Iteration',ylab='Acceptance rate',type='l',
                 main=paste('Acceptance rate:',k,'branches'))
            for(numi in 2:numchain){
                points(sampchaink.accept.rate[[numi]],ylim=c(0,1),xlim=c(1,max(xmax)),col=numi,type='l')
            }
            par(mfrow=c(1,1))
        }
        sampchain[[ki]] = sampchaink

        ki = ki + 1
    }
    if(plot.likelihood) {
        dev.off()
    }
    return(sampchain)
}

###################### canopy.cluster #####################################
canopy.cluster.bugfix <- function(R, X, num_cluster, num_run, Mu.init = NULL,
                           Tau_Kplus1 = NULL) {
  if (is.null(Tau_Kplus1)) {
    Tau_Kplus1 = 0
  }
  VAF = R/X
  s = nrow(R)
  r = pmax(R, 1)
  x = pmax(X, 1)
  Mu_output = Tau_output = pGrank_output = bic_output = vector("list",
                                                               length(num_cluster))
  for (K in num_cluster) {
    cat("Running EM with", K, "clusters...\t")
    Mu_run = Tau_run = pGrank_run = bic_run = vector("list",
                                                     num_run)
    for (run in 1:num_run) {
      cat(run, "  ")
      bic.temp = 0
      Tau = rep(NA, K + 1)
      Tau[K + 1] = Tau_Kplus1
      Tau[1:K] = (1 - Tau_Kplus1)/K
      if (K == 1) {
        Mu = t(as.matrix(apply(R/X, 2, mean)))
      }
      else {
        if (run == 1 & (!is.null(Mu.init))) {
          Mu = Mu.init
        }
        else if (run <= (num_run/2)) {
          VAF.pheat = pheatmap(VAF, cluster_rows = TRUE,
                               cluster_cols = FALSE, kmeans_k = K, silent = TRUE,
                               clustering_distance_rows = "euclidean")
          Mu = pmax(VAF.pheat$kmeans$centers, 0.001)
        }
        else {
          if (ncol(R) > 1) {
            VAF.pheat = pheatmap(VAF, cluster_rows = TRUE,
                                 cluster_cols = FALSE, kmeans_k = K, silent = TRUE,
                                 clustering_distance_rows = "correlation")
            Mu = pmax(VAF.pheat$kmeans$centers, 0.001)
          }
          else {
            VAF.pheat = pheatmap(VAF, cluster_rows = TRUE,
                                 cluster_cols = FALSE, kmeans_k = K, silent = TRUE,
                                 clustering_distance_rows = "euclidean")
            Mu = pmax(VAF.pheat$kmeans$centers, 0.001)
          }
        }
      }
      diff = 1
      numiters = 1
      while (diff > 0.001 || numiters <= 30) {
        numiters = numiters + 1
        pG = canopy.cluster.Estep.bugfix(Tau, Mu, r, x)
        curM = canopy.cluster.Mstep.bugfix(pG, R, X, Tau_Kplus1)
        curTau = curM$Tau
        curMu = curM$Mu
        diff = max(max(abs(Tau - curTau)), max(abs(Mu -
                                                     curMu)))
        Mu = curMu
        Tau = curTau
      }
      dim(pG)
      pGrank = apply(pG, 2, which.max)
      for (i in 1:s) {
        if (pGrank[i] <= K) {
          muk = Mu[pGrank[i], ]
          for (j in 1:ncol(R)) {
            bic.temp = bic.temp + log(Tau[pGrank[i]]) +
              r[i, j] * log(muk[j]) + (x[i, j] - r[i,
                                                   j]) * log(1 - muk[j])
          }
        }
        if (pGrank[i] == (K + 1)) {
          for (j in 1:ncol(R)) {
            bic.temp = bic.temp + log(Tau[pGrank[i]]) +
              lbeta(r[i, j] + 1, x[i, j] - r[i, j] +
                      1)
          }
        }
      }
      bic.temp = 2 * bic.temp - 3 * (length(Tau) - 2 +
                                       length(Mu)) * log(length(R) + length(X))
      Mu_run[[run]] = Mu
      Tau_run[[run]] = Tau
      pGrank_run[[run]] = pGrank
      bic_run[[run]] = bic.temp
    }
    Mu_output[[which(num_cluster == K)]] = Mu_run[[which.max(bic_run)]]
    Tau_output[[which(num_cluster == K)]] = Tau_run[[which.max(bic_run)]]
    pGrank_output[[which(num_cluster == K)]] = pGrank_run[[which.max(bic_run)]]
    bic_output[[which(num_cluster == K)]] = bic_run[[which.max(bic_run)]]
    cat("\n")
  }
  bic_output = as.numeric(bic_output)
  Mu = round(Mu_output[[which.max(bic_output)]], 3)
  Tau = round(Tau_output[[which.max(bic_output)]], 3)
  pGrank = pGrank_output[[which.max(bic_output)]]
  sna_cluster = pGrank
  return(list(bic_output = bic_output, Mu = Mu, Tau = Tau,
              sna_cluster = sna_cluster))
}

canopy.cluster.Mstep.bugfix <- function(pG, R, X, Tau_Kplus1) {
  s = nrow(R)
  K = nrow(pG) - 1
  Tau = rep(NA, K + 1)
  Tau[1:K] = (1 - Tau_Kplus1) * apply(pG[1:K,, drop=FALSE], 1, sum)/(s -
                                                            sum(pG[K + 1, ]))
  Tau[K + 1] = Tau_Kplus1
  pGtemp = pG[1:K,, drop=FALSE]
  Mu = (pGtemp %*% R)/(pGtemp %*% X)
  Mu = round(pmax(Mu, 1e-04), 4)
  return(list(Mu = Mu, Tau = Tau))
}

canopy.cluster.Estep.bugfix <- function(Tau, Mu, R, X) {
  s = nrow(R)
  K = nrow(Mu)
  Mu = pmax(Mu, 0.001)
  pG = matrix(nrow = K + 1, ncol = s, data = log(Tau))
  pG[1:K,] = pG[1:K,] + log(Mu) %*% t(R) + log(1 - Mu) %*%
    (t(X - R))
  for (j in 1:ncol(R)) {
    pG[K + 1, ] = pG[K + 1, ] + lbeta(R[, j] + 1, X[, j] -
                                        R[, j] + 1)
  }
  if (Tau[length(Tau)] != 0) {
    pGtemp = pG
    pGtemp = pGtemp - matrix(ncol = ncol(pGtemp), nrow = nrow(pGtemp),
                             data = apply(pGtemp, 2, max), byrow = TRUE)
    pGtemp = exp(pGtemp)
    pGtemp = pGtemp/matrix(nrow = nrow(pG), ncol = ncol(pG),
                           data = colSums(pGtemp), byrow = TRUE)
    pG[K + 1, (rank(pGtemp[K + 1, ], ties.method = "random")) <=
         (ncol(pG) * (1 - Tau[K + 1]))] = -Inf
  }
  pG = pG - matrix(ncol = ncol(pG), nrow = nrow(pG), data = apply(pG,
                                                                  2, max), byrow = TRUE)
  pG = exp(pG)
  pG = pG/matrix(nrow = nrow(pG), ncol = ncol(pG), data = colSums(pG),
                 byrow = TRUE)
  return(pG)
}
