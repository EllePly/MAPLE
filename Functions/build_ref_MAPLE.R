# Reference network construction functions for MAPLE
#
# These functions build protein-pair reference adjacency matrices from
# protein complex annotations such as CORUM or user-defined complex lists.


# Convert a wide-format complex table into a named complex list.
#
# Input format:
#   Each column is one complex.
#   Each non-empty value within a column is one protein/member.
wide_complex_table_to_list_maple <- function(
    complex_table
) {
  
  complex_list <- lapply(complex_table, function(proteins) {
    proteins <- proteins[!is.na(proteins) & proteins != ""]
    unique(as.character(proteins))
  })
  
  names(complex_list) <- colnames(complex_table)
  
  complex_list <- complex_list[lengths(complex_list) > 0]
  
  return(complex_list)
}


# Map protein identifiers in a complex list.
#
# Example:
#   UniProt IDs -> Gene names
map_complex_ids_maple <- function(
    complex_list,
    id_map,
    from_col = "UniprotID",
    to_col = "GeneName",
    keep_unmapped = TRUE
) {
  
  if (!from_col %in% colnames(id_map)) {
    stop("from_col not found in id_map: ", from_col)
  }
  
  if (!to_col %in% colnames(id_map)) {
    stop("to_col not found in id_map: ", to_col)
  }
  
  id_to_name <- setNames(
    id_map[[to_col]],
    id_map[[from_col]]
  )
  
  mapped_list <- lapply(complex_list, function(ids) {
    
    ids <- as.character(ids)
    mapped <- id_to_name[ids]
    
    if (isTRUE(keep_unmapped)) {
      mapped <- ifelse(is.na(mapped), ids, mapped)
    } else {
      mapped <- mapped[!is.na(mapped)]
    }
    
    mapped <- unique(unname(mapped))
    mapped <- mapped[!is.na(mapped) & mapped != ""]
    
    return(mapped)
  })
  
  mapped_list <- mapped_list[lengths(mapped_list) > 0]
  
  return(mapped_list)
}


# Build a symmetric protein-pair adjacency matrix from a named complex list.
#
# Proteins belonging to the same complex are labelled as 1.
# All other protein pairs are labelled as 0.
build_reference_adjacency_from_complex_list_maple <- function(
    complex_list,
    remove_empty = TRUE,
    remove_singletons = TRUE
) {
  
  if (isTRUE(remove_empty)) {
    complex_list <- lapply(complex_list, function(x) {
      x <- unique(as.character(x))
      x <- x[!is.na(x) & x != ""]
      x
    })
    
    complex_list <- complex_list[lengths(complex_list) > 0]
  }
  
  if (isTRUE(remove_singletons)) {
    complex_list <- complex_list[lengths(complex_list) > 1]
  }
  
  proteins <- sort(unique(unlist(complex_list)))
  
  if (length(proteins) == 0) {
    stop("No proteins found in complex_list.")
  }
  
  ref_mat <- matrix(
    0,
    nrow = length(proteins),
    ncol = length(proteins),
    dimnames = list(proteins, proteins)
  )
  
  for (complex_name in names(complex_list)) {
    
    members <- intersect(complex_list[[complex_name]], proteins)
    
    if (length(members) < 2) {
      next
    }
    
    ref_mat[members, members] <- 1
  }
  
  diag(ref_mat) <- 0
  
  return(ref_mat)
}