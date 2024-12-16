#!/usr/bin/env Rscript

# combine all coverage.tsv from all samples
# generate an tsv and html table 

# arg[1] is a pattern to list files
# arg[2] is workflow id
library(vroom)
library(dplyr)
library(DT)

arg <- commandArgs(trailingOnly = T)
tsvfiles <- list.files(pattern = arg[1], full.names = T)

cols <- c(
    "qname", "rname", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq" 
)
df <- vroom(file = tsvfiles, col_names = cols, col_select = -c('startpos', 'endpos'))

rowCallback <- c(
  "function(row, data){",
  "  for(var i=0; i<data.length; i++){",
  "    if(data[i] === null){",
  "      $('td:eq('+i+')', row).html('-')",
  "        .css({'color': 'rgb(151,151,151)', 'font-style': 'italic'});",
  "    }",
  "  }",
  "}"  
)
# lexocographical arrange of file
locale <- list(locale = "en_US", numeric = TRUE)

finaltable <- 
  DT::datatable(
    dplyr::arrange(df, stringi::stri_rank(qname, opts_collator = locale)),
    class = 'compact',
    # caption = paste0("Run name: ", arg[2], " | Time: ", format.POSIXct(Sys.time())),
    caption = htmltools::tags$caption(
      style = 'caption-side: bottom; text-align: left; color: grey;',
      paste0(arg[2], " | ", format.POSIXct(Sys.time()))
    ),
    # style = 'bootstrap',
    escape = F, filter = 'top',
    extensions = 'Buttons', rownames = FALSE,
    options = list(
      searchHighlight = TRUE,
      rowCallback = JS(rowCallback),
      autoWidth = TRUE, pageLength = 125,
      dom = 'Btp',
      paging = FALSE,
      buttons = c('copy', 'csv', 'excel')
    )
  ) %>% 
  DT::formatRound('coverage', digits = 2) %>%
  DT::formatStyle('coverage', color = styleInterval(c(50, 80), c('#e74c3c', '#f5b041', '#1e8449'))) %>%
  DT::formatStyle('meanbaseq', color = styleInterval(c(20, 25), c('#e74c3c', '#f5b041', '#1e8449'))) %>%
  DT::formatStyle('meanmapq', color = styleInterval(c(40, 50), c('#e74c3c', '#f5b041', '#1e8449')))

#write.csv(df, file = '00-alignment-summary.tsv', sep = '\t', row.names = F, col.names = T)
DT::saveWidget(finaltable, '00-alignment-summary.html')