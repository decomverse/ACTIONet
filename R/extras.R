
map.clusters <- function(Labels, clusters) {
    N = length(Labels)
    cluster.ids = sort(unique(clusters))
    if (is.factor(Labels)) 
        Label.ids = levels(Labels) else Label.ids = sort(unique(Labels))
    
    W = sapply(Label.ids, function(label) {
        idx1 = which(Labels == label)
        n1 = length(idx1)
        if (n1 == 0) 
            return(array(0, length(cluster.ids)))
        
        log.pvals = sapply(cluster.ids, function(cluster.id) {
            idx2 = which(clusters == cluster.id)
            n2 = length(idx2)
            if (n2 == 0) 
                return(0)
            
            success = intersect(idx1, idx2)
            
            pval = phyper(length(success) - 1, n1, N - n1, n2, lower.tail = FALSE)
            return(-log10(pval))
        })
        return(log.pvals)
    })
    
    W[is.na(W)] = 0
    W[W > 300] = 300
    
    W.matched = MWM_hungarian(W)
    W.matched = as(W.matched, "dgTMatrix")
    
    updated.Labels = rep(NA, length(Labels))
    matched.clusters = W.matched@i + 1
    matched.celltype = Label.ids[W.matched@j + 1]
    for (k in 1:length(matched.clusters)) {
        updated.Labels[clusters == matched.clusters[k]] = matched.celltype[[k]]
    }
    
    return(updated.Labels)
}


orthoProject <- function(A, S) {
    A = scale(A)
    S = scale(S)
    A_r = A - S %*% MASS::ginv(t(S) %*% S) %*% (t(S) %*% A)
    A_r = scale(A_r)
    return(A_r)
}

is.sparseMatrix <- function(aa) {
    return(length(which(is(aa) == "sparseMatrix")) != 0)
}

# HGT tail bound
Kappa <- function(p, q) {
    kl = array(1, length(p))
    
    suppressWarnings({
        a = p * log(p/q)
    })
    a[p == 0] = 0
    
    suppressWarnings({
        b = (1 - p) * log((1 - p)/(1 - q))
    })
    b[p == 1] = 0
    
    k = a + b
    return(k)
}

HGT_tail <- function(population.size, success.count, sample.size, observed.success) {
    if (sum(success.count) == 0) 
        return(rep(0, success.count))
    
    success.rate = success.count/population.size
    expected.success = sample.size * success.rate
    delta = (observed.success/expected.success) - 1
    
    log.tail_bound = sample.size * Kappa((1 + delta) * success.rate, success.rate)
    log.tail_bound[delta < 0] = 0
    log.tail_bound[is.na(log.tail_bound)] = 0
    
    return(log.tail_bound)
}

# define utility function to adjust fill-opacity using css
fillOpacity <- function(., alpha = 0.5) {
    css <- sprintf("<style> .js-fill { fill-opacity: %s !important; } </style>", alpha)
    prependContent(., HTML(css))
}

mycircle <- function(coords, v = NULL, params) {
    library(igraph)
    
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
    }
    vertex.size <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
    }
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
        vertex.frame.color <- vertex.frame.color[v]
    }
    vertex.frame.width <- params("vertex", "frame.width")
    if (length(vertex.frame.width) != 1 && !is.null(v)) {
        vertex.frame.width <- vertex.frame.width[v]
    }
    
    mapply(coords[, 1], coords[, 2], vertex.color, vertex.frame.color, vertex.size, vertex.frame.width, FUN = function(x, y, bg, fg, size, 
        lwd) {
        symbols(x = x, y = y, bg = bg, fg = fg, lwd = lwd, circles = size, add = TRUE, inches = FALSE)
    })
}

DECODE = function(hash_str, settings) {
    if (hash_str == "") 
        stop("decode: invalid hashid")
    
    salt = settings$salt
    alphabet = settings$alphabet
    separator = settings$separator
    guards = settings$guards
    
    parts = SPLIT(hash_str, guards)
    hashid = ifelse(2 <= length(parts) & length(parts) <= 3, parts[2], parts[1])
    
    if (hashid == "") 
        stop("decode: invalid hashid, cannot decode")
    lottery = substr(hashid, 1, 1)
    hashid = substr(hashid, 2, nchar(hashid))
    
    hash_parts = SPLIT(hashid, separator)
    unhashed_parts = c()
    for (p in hash_parts) {
        alphabet_salt = substr(paste0(lottery, salt, alphabet), 1, nchar(alphabet))
        alphabet = shuffle(alphabet, alphabet_salt)
        unhashed_parts = c(unhashed_parts, unhash(p, alphabet))
    }
    
    rehash = tryCatch({
        ENCODE(unhashed_parts, settings)
    }, error = function(e) {
        stop("decode: invalid hashid, cannot decode")
    })
    if (!all(hash_str == rehash)) {
        stop("decode: invalid hashid, cannot decode")
    }
    
    return(unhashed_parts)
}

ENCODE = function(int, settings) {
    if (!all(c("alphabet", "salt", "guards", "separator", "min_length") %in% names(settings))) {
        stop("encode: missing some parameters in settings list")
    }
    if (any(int < 0)) {
        stop("encode: numbers must be non-negative")
    }
    if (any(int%%1 != 0)) {
        stop("encode: numbers must be integers")
    }
    if (length(int) < 1) {
        stop("encode: Invalid length!")
    }
    
    alphabet = settings$alphabet
    salt = settings$salt
    guards = settings$guards
    separator = settings$separator
    min_length = settings$min_length
    alphabet_len = nchar(settings$alphabet)
    sep_len = nchar(settings$separator)
    
    vec_hash = sum(sapply(1:length(int), function(i) {
        int[i]%%(100 + i - 1)
    }))
    # lottery character
    lottery = substr(alphabet, (vec_hash%%alphabet_len) + 1, (vec_hash%%alphabet_len) + 1)
    encoded = lottery
    
    for (i in 1:length(int)) {
        alphabet_salt = substr(paste0(lottery, salt, alphabet), 1, alphabet_len)
        alphabet = shuffle(alphabet, alphabet_salt)
        last = hash(int[i], alphabet)
        encoded = paste0(encoded, last)
        int[i] = int[i]%%(ascii_val(substr(last, 1, 1)) + (i - 1))
        encoded = paste0(encoded, substr(separator, (int[i]%%sep_len + 1), (int[i]%%sep_len + 1)))
    }
    
    encoded = substr(encoded, 1, nchar(encoded) - 1)
    if (nchar(encoded) <= min_length) {
        encoded = enforce_min_length(encoded, min_length, alphabet, guards, vec_hash)
    }
    
    return(encoded)
}

hashid_settings = function(salt, min_length = 0, alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890", sep = "cfhistuCFHISTU") {
    
    alphabet_vec = unique(strsplit(alphabet, split = "")[[1]])
    sep_vec = unique(strsplit(sep, split = "")[[1]])
    
    separator_ = paste(intersect(sep_vec, alphabet_vec), collapse = "")
    alphabet_ = paste(setdiff(alphabet_vec, sep_vec), collapse = "")
    
    if (nchar(separator_) + nchar(alphabet_) < 16) {
        # if(nchar(alphabet_) < 16) {
        stop("hashid_settings: Alphabet must be at least 16 unique characters.")
    }
    
    separator_ = shuffle(separator_, salt)
    min_separators = ceiling(nchar(alphabet_)/3.5)
    
    ## if needed get more separators from alphabet ##
    if (nchar(separator_) < min_separators) {
        if (min_separators == 1) 
            min_separators = 2
        split_at = min_separators - nchar(separator_)
        separator_ = paste0(separator_, substr(alphabet_, 1, split_at))
        alphabet_ = substr(alphabet_, split_at + 1, nchar(alphabet_))
    }
    
    alphabet_ = shuffle(alphabet_, salt)
    num_guards = ceiling(nchar(alphabet_)/12)
    
    if (nchar(alphabet_) < 3) {
        guards_ = substring(separator_, 1, num_guards)
        separator_ = substr(separator_, num_guards + 1, nchar(separator_))
    } else {
        guards_ = substring(alphabet_, 1, num_guards)
        alphabet_ = substr(alphabet_, num_guards + 1, nchar(alphabet_))
    }
    
    return(list(alphabet = alphabet_, salt = salt, guards = guards_, separator = separator_, min_length = min_length))
}

ascii_val = function(char) {
    if (!is.character(char)) 
        stop("ascii_val: must be character")
    strtoi(charToRaw(char), 16)
}

base16_to_dec = function(str_16) {
    str_vec = strsplit(tolower(str_16), split = "")[[1]]
    str_vec = sapply(str_vec, function(x) {
        if (x %in% as.character(0:9)) {
            as.numeric(x)
        } else if (x %in% c("a", "b", "c", "d", "e", "f")) {
            ascii_val(x) - 87
        } else {
            stop("base16_to_dec: Invalid hex character")
        }
    })
    
    vec_pwrs = 16^(rev(1:length(str_vec)) - 1)
    
    sum(vec_pwrs * str_vec)
}

dec_to_base16 = function(dec) {
    num_vec = c()
    while (dec > 0) {
        rem = dec%%16
        num_vec = c(rem, num_vec)
        dec = floor(dec/16)
    }
    
    hex_vec = sapply(num_vec, function(x) {
        if (x < 10) {
            return(x)
        } else {
            base::letters[x - 9]
        }
    })
    
    paste(hex_vec, collapse = "")
}

enforce_min_length = function(encoded, min_length, alphabet, guards, values_hash) {
    
    guards_len = nchar(guards)
    guards_idx = (values_hash + ascii_val(substr(encoded, 1, 1)))%%guards_len + 1
    encoded = paste0(substr(guards, guards_idx, guards_idx), encoded)
    
    if (nchar(encoded) < min_length) {
        guards_idx = (values_hash + ascii_val(substr(encoded, 3, 3)))%%guards_len + 1
        encoded = paste0(encoded, substr(guards, guards_idx, guards_idx))
    }
    
    split_at = nchar(alphabet)/2 + 1
    while (nchar(encoded) < min_length) {
        alphabet = shuffle(alphabet, alphabet)
        encoded = paste0(substr(alphabet, split_at, nchar(alphabet)), encoded, substr(alphabet, 1, split_at - 1))
        excess = nchar(encoded) - min_length
        
        if (excess > 0) {
            from_index = floor(excess/2) + 1
            encoded = substr(encoded, from_index, from_index + min_length - 1)
        }
    }
    
    return(encoded)
}

hash = function(number, alphabet) {
    alphabet_len = nchar(alphabet)
    alphabet_vec = strsplit(alphabet, split = "")[[1]]
    
    hashed = c()
    while (number > 0) {
        hash_idx = (number%%alphabet_len) + 1
        hashed = c(alphabet_vec[hash_idx], hashed)
        number = floor(number/alphabet_len)
    }
    
    return(paste(hashed, collapse = ""))
}

shuffle = function(string, salt) {
    salt_len = nchar(salt)
    str_len = nchar(string)
    if (salt_len < 1 | str_len < 2) 
        return(string)
    
    salt_sum = 0
    salt_index = 1
    string_vec = strsplit(string, split = "")[[1]]
    salt_vec = strsplit(salt, split = "")[[1]]
    for (i_str in rev(2:str_len)) {
        ## Pseudo Randomize based on salt ##
        salt_index = (salt_index - 1)%%(salt_len) + 1
        salt_int_val = ascii_val(salt_vec[salt_index])
        salt_sum = salt_sum + salt_int_val
        swap_pos = (salt_sum + salt_index + salt_int_val - 1)%%(i_str - 1) + 1
        
        ## Swap positions ##
        temp = string_vec[swap_pos]
        string_vec[swap_pos] = string_vec[i_str]
        string_vec[i_str] = temp
        salt_index = salt_index + 1
    }
    
    return(paste(string_vec, collapse = ""))
}

SPLIT = function(string, splitters) {
    string_vec = strsplit(string, split = "")[[1]]
    split_vec = strsplit(splitters, split = "")[[1]]
    
    word = ""
    words = c()
    for (i in 1:length(string_vec)) {
        if (string_vec[i] %in% split_vec) {
            words = c(words, word)
            word = ""
        } else {
            word = paste0(word, string_vec[i])
        }
    }
    words = c(words, word)
    
    return(words)
}

unhash = function(hashed, alphabet) {
    hashed_len = nchar(hashed)
    alphabet_len = nchar(alphabet)
    alphabet_vec = strsplit(alphabet, split = "")[[1]]
    hashed_vec = strsplit(hashed, split = "")[[1]]
    
    number = 0
    for (i in 1:hashed_len) {
        position = which(alphabet_vec == hashed_vec[i]) - 1
        number = number + (position * alphabet_len^(hashed_len - i))
    }
    
    return(number)
}



combine.logPvals <- function(logPvals, top.len = NULL, base = 10) {
    if (is.null(top.len)) {
        top.len = nrow(logPvals)
    }
    kappa = 1/log(exp(1), base = base)
    logPvals = kappa * logPvals
    
    combbined.log.pvals = -apply(logPvals, 2, function(lx) {
        perm = order(lx, decreasing = T)
        
        return(log(top.len) - logSumExp(lx[perm[1:top.len]]))
    })
    
    return(combbined.log.pvals)
}



reannotate.labels <- function(ace, Labels) {
    Labels = preprocess.labels(Labels, ace)
    
    
    Annot = sort(unique(Labels))
    idx = match(Annot, Labels)
    names(Annot) = names(Labels)[idx]
    
    new.Labels = match(names(Labels), names(Annot))
    names(new.Labels) = as.character(names(Labels))
    
    return(new.Labels)
}

gen.colors <- function(Pal, color.no, plot.cols) {
    color.no = min(color.no, length(Pal) - 1)
    colors.RGB = t(col2rgb(Pal)/256)
    colors.Lab = grDevices::convertColor(color = colors.RGB, from = "sRGB", to = "Lab")
    
    set.seed(0)
    W0 = t(kmeans(colors.Lab, color.no)$centers)
    AA.out = ACTIONet::runAA(X, W0)
    
    C = AA.out$C
    arch.colors = t(t(colors.Lab) %*% C)
    
    new.colors = rgb(grDevices::convertColor(color = arch.colors, from = "Lab", to = "sRGB"))
    if (plot.cols) 
        scales::show_col(new.colors)
    return(new.colors)
}


doubleNorm <- function(Enrichment, log.transform = T, min.threshold = 0) {
    # if(min(Enrichment) < 0) { Enrichment = exp(Enrichment) }
    Enrichment[Enrichment < min.threshold] = 0
    
    if ((max(Enrichment) > 100) & (log.transform == T)) {
        Enrichment = log1p(Enrichment)
    }
    Enrichment[is.na(Enrichment)] = 0
    
    rs = sqrt(Matrix::rowSums(Enrichment))
    rs[rs == 0] = 1
    D_r = Matrix::Diagonal(nrow(Enrichment), 1/rs)
    
    cs = sqrt(Matrix::colSums(Enrichment))
    cs[cs == 0] = 1
    D_c = Matrix::Diagonal(ncol(Enrichment), 1/cs)
    
    Enrichment.scaled = as.matrix(D_r %*% Enrichment %*% D_c)
    
    Enrichment.scaled = Enrichment.scaled/max(Enrichment.scaled)
    return(Enrichment.scaled)
}

assess.label.local.enrichment <- function(P, Labels) {
    if (is.null(names(Labels))) {
        names(Labels) = as.character(Labels)
    }
    counts = table(Labels)
    p = counts/sum(counts)
    Annot = names(Labels)[match(as.numeric(names(counts)), Labels)]
    
    X = sapply(names(p), function(label) {
        x = as.numeric(Matrix::sparseVector(x = 1, i = which(Labels == label), length = length(Labels)))
    })
    colnames(X) = Annot
    
    Exp = array(1, nrow(P)) %*% t(p)
    Obs = as(P %*% X, "dgTMatrix")
    
    # Need to rescale due to missing values within the neighborhood
    rs = Matrix::rowSums(Obs)
    Obs = Matrix::sparseMatrix(i = Obs@i + 1, j = Obs@j + 1, x = Obs@x/rs[Obs@i + 1], dims = dim(Obs))
    
    Lambda = Obs - Exp
    
    
    w2 = Matrix::rowSums(P^2)
    Nu = w2 %*% t(p)
    
    a = as.numeric(qlcMatrix::rowMax(P)) %*% t(array(1, length(p)))
    
    
    logPval = (Lambda^2)/(2 * (Nu + (a * Lambda)/3))
    logPval[Lambda < 0] = 0
    logPval[is.na(logPval)] = 0
    
    logPval = as.matrix(logPval)
    
    colnames(logPval) = Annot
    
    
    max.idx = apply(logPval, 1, which.max)
    updated.Labels = as.numeric(names(p))[max.idx]
    names(updated.Labels) = Annot[max.idx]
    
    updated.Labels.conf = apply(logPval, 1, max)
    
    res = list(Labels = updated.Labels, Labels.confidence = updated.Labels.conf, Enrichment = logPval)
    
    return(res)
}

make.fname <- function(str) {
    nn = stringr::str_replace(str, " ", "_")
    nn = stringr::str_replace(nn, fixed("."), "_")
    nn = stringr::str_replace(nn, fixed("/"), "_")
    
    return(nn)
}

add.count.metadata <- function(ace) {
    ace$n_counts = Matrix::colSums(counts(ace))
    ace$n_genes = Matrix::colSums(counts(ace) > 0)
    rowData(ace)$n_cells = Matrix::rowSums(counts(ace) > 0)
    
    return(ace)
}
