# VCF validator based on 
# https://samtools.github.io/hts-specs/VCFv4.2.pdf


# helper functions #####################
log <-function(i, text){
  paste("line ",i,": ",text, sep = "")
}
isFloat <- function(text){
  grepl("^[0-9]+(\\.[0-9]+)?$", text, perl=T)
}
containsNucleotidesRef <- function(text){
  grepl("^([ACGTN]+|\\.$)|(^[acgtn]+|\\.$)", perl=TRUE, text)
}
containsNucleotidesAlt <- function(text){
  grepl("^([ACGTN]+|\\*|<CNV>|\\.)(,([acgtn]+|\\*|<CNV>|\\.))*$", perl = TRUE, text)
}
isInfo <- function(text){
  grepl("^.+(=.+(,.+)*)?(;.+(=.+(,.+)*)?)*$", text)
}
isGenotype <- function(text){
  grepl("^([\\d|\\.]\\|[\\d|\\.])|([\\d|\\.]/[\\d|\\.])$", perl = TRUE, text)
}

parseInfo <- function(text){
  keys = c()
  tokens = unlist(strsplit(text, ";"))
  
  for(i in 1:length(tokens)){
    keys <- c(keys, unlist(strsplit(tokens[i][1],'='))[1])
  }
  return(keys)
}
#######################################

.validate_vcf <- function(vcf_file, check_genotype=FALSE, AF_softcheck=TRUE){
  ## Read vcf in memory
  col_num = utils::read.table(file = vcf_file, header = F, stringsAsFactors = FALSE,
                      blank.lines.skip=T, comment.char = "#", fill = T,
                      sep="\t",nrows = 1)
  
  f <- utils::read.table(file = vcf_file, header = F, stringsAsFactors = FALSE,
                  blank.lines.skip=T, comment.char = "", fill = T,
                  sep="\t",col.names = paste("V",1:length(col_num)))
  rm(col_num)
  
  ## Output data structures
  warnings =c()
  errors =c()
  
  ## First check: File format
  if(!startsWith(f[1,1], "##fileformat")){
    errors = c(errors, log(1, "The vcf file must start with ##fileformat"))
  }
  
  ## Load metadata
  info_id = c()
  filter_id = c("PASS")
  format_id = c()
  samples_num = 0
  variants_list = c()
  AF_in_FORMAT = FALSE
  AF_format_max_warnings = 1
  AF_in_INFO = FALSE
  AF_info_max_warnings = 1
  AF_info_format_max_warnings = 1
  AF_info_format_both_max_warnings = 1
  AF_info_format_either_max_warnings = 1
  REF_ERROR = 1
  ALT_ERROR = 1
  REF_ALT_max_errors = 1
  ALT_REF_max_errors = 1
  
  for(i in 1:nrow(f)){
    if(startsWith(f[i,1], "##")){
      if(startsWith(f[i,1],"##INFO")){
        # INFO fields should be described as: first four keys are required, source and version are recommended
        # ##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">
        temptokens = unlist(strsplit(gsub("##INFO=<","",f[i,1]),","))
        if(length(temptokens)<4){
          errors = c(errors, log(i, "The INFO field must contain at least the 'ID','Number','Type','Description' keys"))
        }else{
          tempkeys = c()
          for(j in 1:length(temptokens)){
            tempsubtokens = unlist(strsplit(temptokens[j], "="))
            tempkeys = c(tempkeys, tempsubtokens[1])
            if(tempsubtokens[1]=="ID"){
              info_id <- c(info_id, tempsubtokens[2])
            }
          }
          if(is.na(match('ID',tempkeys))){  #ID
            errors = c(errors, log(i, "The INFO field must contain the 'ID' key"))
          }
          if(is.na(match('Number',tempkeys))){  #Number
            errors = c(errors, log(i, "The INFO field must contain the 'Number' key"))
          }
          if(is.na(match('Type',tempkeys))){  #Type
            errors = c(errors, log(i, "The INFO field must contain the 'Type' key"))
          }
          if(is.na(match('Description',tempkeys))){  #Description
            errors = c(errors, log(i, "The INFO field must contain the 'Description' key"))
          }
          # check if AF is present throughout the INFO header lines, if requested
          if(AF_softcheck==TRUE){
            if('AF' %in% info_id){
              AF_in_INFO = TRUE  # initialized to FALSE
            }
          }
          # ----------------------------------------------------------
        }
  
      } 
      else if(startsWith(f[i,1],"##FILTER")){
        # FILTERs that have been applied to the data should be described as follows:
        # ##FILTER=<ID=ID,Description="description">
        temptokens = unlist(strsplit(gsub("##FILTER=<","",f[i,1]),","))
        tempkeys = c()
        for(j in 1:length(temptokens)){
          tempsubtokens = unlist(strsplit(temptokens[j], "="))
          tempkeys = c(tempkeys, tempsubtokens[1])
          if(tempsubtokens[1]=="ID"){
            filter_id <- c(filter_id, tempsubtokens[2])
          }
        }
        if(is.na(match('ID',tempkeys))){  #ID
          errors = c(errors, log(i, "The FILTER field must contain the 'ID' key"))
        }
        if(is.na(match('Description',tempkeys))){  #Description
          errors = c(errors, log(i, "The FILTER field must contain the 'Description' key"))
        }
      } 
      else if(startsWith(f[i,1],"##FORMAT")){
        # Genotype fields specified in the FORMAT field should be described as follows:
        # ##FORMAT=<ID=ID,Number=number,Type=type,Description="description">
        temptokens = unlist(strsplit(gsub("##FORMAT=<","",f[i,1]),","))
        if(length(temptokens)<4){
          errors = c(errors, log(i, "The FORMAT field must contain at least the 'ID','Number','Type','Description' keys"))
        }else{
          tempkeys = c()
          for(j in 1:length(temptokens)){
            tempsubtokens = unlist(strsplit(temptokens[j], "="))
            tempkeys = c(tempkeys, tempsubtokens[1])
            if(tempsubtokens[1]=="ID"){
              format_id <- c(format_id, tempsubtokens[2])
            }
          }
          if(is.na(match('ID',tempkeys))){  #ID
            errors = c(errors, log(i, "The FORMAT field must contain the 'ID' key"))
          }
          if(is.na(match('Number',tempkeys))){  #Number
            errors = c(errors, log(i, "The FORMAT field must contain the 'Number' key"))
          }
          if(is.na(match('Type',tempkeys))){  #Type
            errors = c(errors, log(i, "The FORMAT field must contain the 'Type' key"))
          }
          if(is.na(match('Description',tempkeys))){  #Description
            errors = c(errors, log(i, "The FORMAT field must contain the 'Description' key"))
          }
          # check if AF is present throughout the FORMAT header lines, if requested
          if(AF_softcheck==TRUE){
            if('AF' %in% format_id){
              AF_in_FORMAT = TRUE  # initialized to FALSE
            }
          }
          # ----------------------------------------------------------
        }
      }
    }
    else if(startsWith(f[i,1], "#")){
      temp_header = tolower(f[i,])
      #Ex. header: CHROM POS ID REF ALT QUAL FILTER INFO <FORMAT NA00001 NA00002 NA00003>
      if(length(format_id)>0){
        samples_num = length(temp_header) - 9
        
        ## First-bis check: presence of multiple samples
        if(samples_num > 1){
          warnings= c(warnings, log(i, paste("The VCF file contains ", 
                                             samples_num, 
                                             " samples. Only the first (", 
                                             temp_header[10],
                                             ") will be used", sep = "")))
        }
      }
      
      ## Second check: presence of required header fields
      chr_idx = match("#chrom", temp_header)
      pos_idx = match("pos", temp_header)
      id_idx = match("id", temp_header)
      ref_idx = match("ref", temp_header)
      alt_idx = match("alt", temp_header)
      qual_idx = match("qual", temp_header)
      filter_idx = match("filter", temp_header)
      info_idx = match("info", temp_header)
      format_idx = match("format", temp_header)
      
      missing_header_item = c()
      if(is.na(chr_idx)){
        missing_header_item = c(missing_header_item, "CHROM")  
      } 
      if(is.na(pos_idx)){
        missing_header_item = c(missing_header_item, "POS")  
      } 
      if(is.na(id_idx)){
        missing_header_item = c(missing_header_item, "ID")  
      } 
      if(is.na(ref_idx)){
        missing_header_item = c(missing_header_item, "REF")  
      } 
      if(is.na(alt_idx)){
        missing_header_item = c(missing_header_item, "ALT")  
      } 
      if(is.na(qual_idx)){
        missing_header_item = c(missing_header_item, "QUAL")  
      } 
      if(is.na(filter_idx)){
        missing_header_item = c(missing_header_item, "FILTER")  
      } 
      if(is.na(info_idx)){
        missing_header_item = c(missing_header_item, "INFO")  
      } 
      if(is.na(format_idx)){
        if(samples_num > 0 & check_genotype==TRUE){
          missing_header_item = c(missing_header_item, "FORMAT")
        }
      }
      
      if(length(missing_header_item)>0){
        errors= c(errors, log(i, paste("missing mandatory fields:", 
                             paste(missing_header_item,collapse=", "))))
      }
      else{
        ## Third check: correct order of header fields
        consecutive_correct_num= which(diff(c(chr_idx,pos_idx,id_idx,ref_idx,alt_idx,
                                              qual_idx,filter_idx,info_idx))==1)
        if(length(consecutive_correct_num) != 7){
          errors= c(errors, log(i, "The correct order of the columns is: CHROM POS ID REF ALT QUAL FILTER INFO <FORMAT SAMPLE1 SAMPLE2 ...>"))
        }
      }
    } else{
      ## Data lines
      AF_in_variant_INFO = FALSE
      AF_in_variant_FORMAT = FALSE
      
      # CHROM must be a String, no white-space permitted, Required
      if(!is.character(f[i,1]) 
         | nchar(gsub("\\s", "", f[i,1])) != nchar(f[i,1])
         | nchar(f[i,1])==0){
        errors= c(errors, log(i, paste("CHROM must be a not-empty String; no white-space permitted")))
      }
      
      # The reference position is a not-empty integer with the 1st base having position 1
      if(nchar(f[i,2])==0){
        errors= c(errors, log(i, paste("POS field is empty")))
      }
      else if(!isFloat(f[i,2])){
        errors= c(errors, log(i, paste("POS must be an integer value; '", f[i,2], "' found", sep="")))
      } else{
        if(as.numeric(f[i,2])<1){
          errors= c(errors, log(i, paste("POS must be an integer with the 1st base having position 1 or greater; '",
                                         f[i,2], "' found", sep="")))  
        }
      }
      
      # Recording chr + pos as chr_pos to account for multi-allelic variants
      variants_list = c(variants_list, paste(f[i,1],"_", f[i,2], sep = ""))
      
      # ID - identifier: Semi-colon separated list of unique identifiers where available.
      # (String, no white-space permitted)
      if(!is.character(f[i,3]) 
         | nchar(gsub("\\s", "", f[i,3])) != nchar(f[i,3])
         | nchar(f[i,3])==0){
        errors= c(errors, log(i, paste("ID must be a not-empty String; no white-space permitted")))
      }
      
      # REF - reference base(s): Each base must be one of A,C,G,T,N (case insensitive)
      ref_provided = TRUE
      if((!is.character(f[i,4]) | !containsNucleotidesRef(f[i,4])) & REF_ERROR>0){
        errors= c(errors, log(i, paste("REF must be one of A,C,G,T,N,. (case insensitive); '", 
                                       f[i,4], "' found", sep="")))
        ref_provided = FALSE
        REF_ERROR = REF_ERROR -1
      }
      
      # ALT - alternate base(s): Comma separated list of alternate non-reference alleles.
      # Options are base Strings made up of the bases A,C,G,T,N,*, (case insensitive)
      # or an angle-bracketed ID String. If there are no alternative alleles, 
      # then the missing value (.) should be used.
      alt_provided = TRUE
      if((!is.character(f[i,5]) | !containsNucleotidesAlt(f[i,5])) & ALT_ERROR>0){
        errors= c(errors, log(i, paste("ALT must be a comma separated list of alternate alleles (A,C,G,T,N,.,*); '", 
                                       f[i,5], "' found", sep="")))
        alt_provided = FALSE
        ALT_ERROR = ALT_ERROR -1
      }
      
      # CHECK presence of either ALT or REF
      if(ref_provided & !alt_provided){
        if(REF_ALT_max_errors > 0 )
          errors= c(errors, log(i, paste("the ALT allele is not provided for this and possibly other variants. If this is a deletion, it should be represented according to the 1000 Genome project as shown in the link: https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40]", sep="")))
          REF_ALT_max_errors = REF_ALT_max_errors -1
      } else if(!ref_provided & alt_provided){
        if(ALT_REF_max_errors > 0 )
          errors= c(errors, log(i, paste("the REF allele is not provided for this and possibly other variants. If this is an insertion, it should be represented according to the 1000 Genome project as shown in the link: https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40]", sep="")))
          ALT_REF_max_errors = ALT_REF_max_errors -1
      }
      
      # QUAL: Phred-scaled quality score for the assertion made in ALT
      if(f[i,6]!="."){
        if(!isFloat(f[i,6])){
          errors= c(errors, log(i, paste("QUAL must be '.' (missing value) or numeric; '", 
                                         f[i,6], "' found", sep="")))
        }
      }
      
      # FILTER - filter status
      # (String, no white-space permitted)
      if(!is.character(f[i,7]) 
         | nchar(gsub("\\s", "", f[i,7])) != nchar(f[i,7])
         | nchar(f[i,7])==0){
        errors= c(errors, log(i, paste("FILTER must be a not-empty String; no white-space permitted")))
      } else{
        temptokens = unlist(strsplit(f[i,7],";"))
        temptokens = temptokens[names(temptokens) != "."] # skip "." key
        if(length(setdiff(temptokens, filter_id)>0)){
          warnings= c(warnings, log(i, paste("Some FILTER elements are not described in the header lines; '", 
                                             paste(setdiff(temptokens, filter_id),collapse=","), "' missing", sep="")))
        }
      } 
      
      # INFO - INFO fields are encoded as a semicolon-separated series of short
      # keys with optional values in the format: <key>=<data>[,data].
      # (String, no white-space permitted)
      # the presence of the AF key is mandatory
      if(!is.character(f[i,8]) 
         | nchar(gsub("\\s", "", f[i,8])) != nchar(f[i,8])
         | nchar(f[i,8])==0){
        errors= c(errors, log(i, paste("INFO must be a not-empty String; no white-space permitted")))
      } else if (!isInfo(f[i,8])){
        errors= c(errors, log(i, paste("INFO field data format not allowed; '", f[i,8], "' found", sep="")))
      } else{
        info_tokens = parseInfo(f[i,8])
        info_tokens = info_tokens[info_tokens != "."] # skip "." key
        missing_info_tokens = setdiff(info_tokens, info_id)
        
        if(length(missing_info_tokens) != 0 ){
          warnings= c(warnings, log(i, paste("Some INFO elements are not described in the header lines; '", 
                                             paste(missing_info_tokens, collapse = ","), "' missing", sep="")))
        }
        # check the presence of AF through the INFO elements, if requested
        if(AF_softcheck==TRUE){
          if(is.na(match('AF',info_tokens)) & AF_in_INFO==TRUE){
            if(AF_info_max_warnings>0){
              warnings= c(warnings, log(i, "AF is described in the header but it is not defined in the INFO field of this and possibly other variants"))
              AF_info_max_warnings = AF_info_max_warnings -1
              AF_in_variant_INFO = FALSE
            }
          } else if(!is.na(match('AF',info_tokens)) & AF_in_INFO==TRUE){
              AF_in_variant_INFO = TRUE
          } else if(is.na(match('AF',info_tokens)) & AF_in_INFO==FALSE){
            AF_in_variant_INFO = FALSE
          } else{
            if(AF_info_max_warnings>0){
              warnings= c(warnings, log(i, "AF is described in the INFO field of this and possibly other variants but not in the the header (INFO)"))
              AF_info_max_warnings = AF_info_max_warnings -1
              AF_in_variant_INFO = TRUE
            }
          }
        }
        # ----------------------------------------------------------------
      }
      
      # FORMAT (if present) - colon-separated alphanumeric String.
      # The first sub-field must always be the genotype (GT) if it is present.
      # There are no required sub-fields.
      if(!is.na(format_idx)){
          subfields = unlist(strsplit(f[i,9], ":"))
          
          if(check_genotype==TRUE & subfields[1] != "GT"){
            errors= c(errors, log(i, paste("The first FORMAT's sub-field must be GT; '", subfields[1], "' found", sep = "")))
          }
          
          # check the presence of AF through the FORMAT elements, if requested
          if(AF_softcheck==TRUE){
            if(is.na(match('AF',subfields)) & AF_in_FORMAT==TRUE){
              if(AF_format_max_warnings>0){
                warnings= c(warnings, log(i, "AF is described in the header but it is not defined in the FORMAT field of this and possibly other variants"))
                AF_format_max_warnings = AF_format_max_warnings-1
                AF_in_variant_FORMAT = FALSE
              }
            }
            else if(!is.na(match('AF',subfields)) & AF_in_FORMAT==TRUE){
              AF_in_variant_FORMAT = TRUE
            }
            else if(is.na(match('AF',subfields)) & AF_in_FORMAT==FALSE){
              AF_in_variant_FORMAT = FALSE
            }
            else{
              if(AF_format_max_warnings>0){
                warnings= c(warnings, log(i, "AF is described in the FORMAT field of this and possibly other variants but not in the the header (FORMAT)"))
                AF_format_max_warnings = AF_format_max_warnings -1
                AF_in_variant_FORMAT = TRUE
              }
            }
          }
          # ----------------------------------------------------------------
          
          # check if all the samples in this vcf have:
          # 1. the same number of subfields than FORMAT
          # 2. their first subfield has a format like number/number or number|number
          # # -> enable this if you want to check the soundness of FORMAT for ALL samples <- for(j in 1:samples_num){
          for(j in 1:1){
            subfields_sample = unlist(strsplit(f[i,9+j], ":"))
            if(length(subfields_sample) != length(subfields)){
              errors= c(errors, log(i, paste("FORMAT and", temp_header[9+j], "must have the same number of sub-fields.", 
                                             temp_header[9+j], "has", length(subfields_sample), "sub-fields;",
                                             "FORMAT has", length(subfields), "sub-fields", sep = " ")))
            }
            
            if(check_genotype==TRUE & !isGenotype(subfields_sample[1])){
              errors= c(errors, log(i, paste(temp_header[9+j], "'s GT subfield must be a legal genotype (e.g., 1/1 or 0|1); '",
                                             subfields_sample[1], "' found", sep = "")))
            }
          } # end for
          
          # check if AF is in INFO and/or FORMAT of this variant, if requested
          if(AF_softcheck==TRUE){
              if(AF_in_variant_INFO==FALSE & AF_in_variant_FORMAT==FALSE & AF_info_format_either_max_warnings>0){
                warnings= c(warnings, log(i, "AF is neither described in FORMAT nor in INFO of this and possibly other variants"))
                AF_info_format_either_max_warnings = AF_info_format_either_max_warnings -1
              } else if(AF_in_variant_INFO==TRUE & AF_in_variant_FORMAT==TRUE & AF_info_format_both_max_warnings>0){
                warnings= c(warnings, log(i, "AF is described both in FORMAT and in INFO of this and possibly other variants. FORMAT's AF will be used"))
                AF_info_format_both_max_warnings = AF_info_format_both_max_warnings -1
              }
          }
      }
    }
  }
  
  # Determine multi-allelic variants
  dup_variants = variants_list[duplicated(variants_list)]
  if(length(dup_variants) > 0){
    for(i in 1:length(dup_variants)){
      var = strsplit(as.character(dup_variants[i]), "_")
      warnings = c(warnings, paste("The multi-allelic variant in chromosome ", var[[1]][1], " and position ", 
                                   var[[1]][2], " cannot be handled by TMBleR", sep=""))
    }
  }


  return(list(warnings=warnings,errors=errors))
}
