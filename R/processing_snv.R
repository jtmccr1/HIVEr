#' Split a file path at /
#'
#' This function splits a string at '/' and returns
#' a one of the entries. You can not return muplitle entries at this point.
#'
#' This is a hepler function in the read_rbind function.
#'
#' @param x The string to be split.
#' @param n The entry to be return
#'
#' @examples
#' split_path("foo/bar/foo.csv",2)


split_path <-function(x,n){
  if(length(n)>1) stop("We can't handel multiple splits now. Sorry")
  split<-strsplit(x,"/")
  split <- unlist(split)
  return(split[n])
}

#' Read in a list of files
#'
#' This function reads in a list of files and combines them into one tibble.
#'
#' There are options to append a run column using a section of the path from
#' each file. When using we should speficy the Id column is a character as some
#' can be interpretted as either a character or a intereger. But others are characters.
#'
#' @param files A list of files to read in and combine
#' @param n the optional entry to append as the run column. It is passed to split_path. It is the
#' nth entry in the path (sparating at '/')
#' @param cols An optional list of column types to pass to read_csv
#'
#'@examples
#'\dontrun{
#'read_rbind(c("../../Analysis/HIVE/data/processed/cali09/all.coverage.csv",
#'"../../Analysis/HIVE/data/processed/cali09_2/all.coverage.csv"),
#'n=7,
#'cols=list(Id=readr::col_character()))
#'}
#'
#' @export
read_rbind<-function(files,n=NULL,cols = NULL){


      data <- dplyr::data_frame(filename = files) %>% # create a data frame
        # holding the file names
        dplyr::mutate(file_contents = purrr::map(filename,~readr::read_csv(.,col_types=cols)) # a new data column
        )
      if(!is.null(n)){
        # add the run name
        data <- data %>% dplyr::mutate(run = purrr::map(filename,~split_path(.,n)))
        data <- data %>% dplyr::select(-filename) %>%
          tidyr::unnest(run) %>%
          tidyr::unnest(file_contents)

      }
      else{
        data <- data %>% dplyr::select(-filename) %>% tidyr::unnest(file_contents)
      }
      return(data)
  }



#' Sift duplicate calls
#' Take in a data frame of one mutation identified in 2 replicates and
#' choose the replicate with the better coverage.
#'
#' @param df a data frame of nrow(2)
#' @examples
#' sift_dups(subset(variants,mutation=="PB1_G530A"))

sift_dups<-function(df){
  if(nrow(df)>2) stop("Too many mutations here")
  df <- dplyr::mutate(df,coverage = cov.tst.bw+cov.tst.fw)
  # this is local coverage from deepsnv
higher_qual<-subset(df,coverage==max(df$coverage))
    if(nrow(higher_qual)>1){ # same cov
      higher_qual<-higher_qual[1,]
    }
    return(higher_qual)
}


#' Identify quality variants
#'
#' This function takes in all variant calls associated with a specimen.
#' If the sample was sequenced in duplicate that must be apparent by a 'run' column
#' that specifies the unique names of the run columns. Samples with gc_ul above 1e5 don't
#' need to be run in duplicate. If the sample has a lower titer but was not sequenced twice this
#' function will through an error.
#'
#' @param df a data frame containing all base calls associated with a specimen or isolate.
#'
#' @examples
#' quality(variants)
#'
#' @export

quality<-function(df){ # this now is written to handle 1 isolate (with possibly 2 sequencing runs at a time)
  #good<-subset(df,gc_ul>=1e5)
  # check only one run/sample
  stopifnot(length(unique(df$LAURING_ID))==1)
  runs <- df %>%
    dplyr::group_by(LAURING_ID) %>%
    dplyr::summarize(runs=length(unique(run)))

  if(runs$runs==2 & unique(df$gc_ul)>1e3){
    if(unique(df$gc_ul)>1e5){
      warning(paste0("Samples ",unique(df$LAURING_ID)), " were sequenced twice even though the titer was ",
              unique(df$gc_ul),".It was  treated as a duplicate in the analysis")
    }

    count_mutations <- df %>% dplyr::group_by(mutation) %>%
      dplyr::summarize(found = length(mutation))

    stopifnot(max(count_mutations$found)<3)

    good_mut<-count_mutations %>% dplyr::filter(found==2)

    good_var <- df %>% dplyr::filter(mutation %in% good_mut$mutation)
# now we only need one instance of the mutation. We'll take the one with the higher frequency.

    dups_good<- good_var %>% dplyr::group_by(mutation)%>%
      dplyr::do(sift_dups(.))

    return(dplyr::ungroup(dups_good))
  }else if(runs$runs==1 & unique(df$gc_ul)>1e5){
    return(df)
  }else{
    warning(paste0(unique(df$LAURING_ID)," Was not sequenced properly given its titer. It will be removed."))
    return(df[F,])
    }
}


#' Correct numerical Ids
#' In our analysis the Id column refers to the sequenced id. In many cases this was
#' a number; but later it matched the SPECID + and extention. Therefore the column
#' should be read as string. However, Pandas for some reason converted the numbers to
#' decimal format. This function corrects that and returns a string number with out the ".0"
#' @param df data frame or tibble to correct
#' @param col column that needs to be corrected.
#'
#' @examples
#' cov_sample<-dplyr::tibble(Id=c("129.0","MH00000"))
#' correct_id(cov_sample,Id)
#' @export

correct_id<-function(df,col){
  quo_col <- enquo(col)
  quo_name <-quo_name(quo_col)

  df<- df %>% dplyr::mutate(!!quo_name := dplyr::if_else(grepl(".0",!!quo_col,fixed=TRUE),
                                           sub(".0","",!!quo_col,fixed = TRUE),
                                           !!quo_col))
}


#' Get sites with diverisity
#'
#' In our analysis we are not interested in loci that only have
#' one allele present within a season. Cutting out such loci
#' helps focus on the important loci and distinctions between samples.
#' This function counts the unique "var" entries in a dataframe according
#' to the set grouping and removes sites that don't meet the cut off.
#'
#' @param df variant data frame
#' @param cutoff integer. The allele must be found more than this many times in the
#' groupings below.
#' @param ... columns to group by. There's no current need for them to be anything but
#' season,chr,pos,pcr_result as this allows us to
#' @return the input data frame but with sites with more alleles than the cut off
#' @examples
#' diverse_sites(small_isnv,1,season,pcr_result,pos,chr)
#'
#' @export
diverse_sites<-function(df,cutoff,...){
  group_var<-rlang::quos(...)
  if(length(unique(df$pcr_result))!=1 |length(unique(df$season))!=1 ){
    warning("You are looking across multiple seasons or strains")
  }
  df_alleles<- df %>% dplyr::group_by(!!!group_var) %>%
    dplyr::summarize(alleles=length(unique(var)))

  df <- dplyr::left_join(df,df_alleles)
  stopifnot(any(!is.na(df$alleles)))
  df <- df %>% dplyr::filter(alleles>cutoff) %>% dplyr::select(-alleles)
}




#' Set monomorphic sites to freq =1
#'
#' In our analysis we calculate frequencies prior to
#' removing false positives. Sometimes this means we have a monomorphic site
#' with the allele present at <100%. I have checked all such sites between 50-90% and
#' in each case this holds true. That work is in 2017_11_13.Rmd- 2017_11_16.Rmd in the
#' notebook directory of the HIVE repo. It is also true we will miss some true variants
#' in the <10% range. Such is life. This script sets all sites that are monomorphic to
#' frequenies of 1.
#'
#' @param df variant data frame
#' @param ... columns to group by. There's no current need for them to be anything but
#' SPECID,season,chr,pos,pcr_result
#' @return the input data frame but with monomorphic frequencies = 1.
#' @examples
#' wacky_isnv<-small_isnv
#' wacky_isnv$freq.var[wacky_isnv$mutation=="PB2_G1602A"]<-0.8
#' monomorphic(wacky_isnv,SPECID,season,chr,pos,pcr_result)
#'
#' @export
monomorphic<-function(df,...){
  group_var<-rlang::quos(...)
  df_alleles<- df %>% dplyr::group_by(!!!group_var) %>%
    dplyr::summarize(alleles=length(unique(var)))

  df <- dplyr::left_join(df,df_alleles)
  stopifnot(any(!is.na(df$alleles)))
  df$freq.var[df$alleles==1]<-1
  df <- df %>% dplyr::select(-alleles)
}
