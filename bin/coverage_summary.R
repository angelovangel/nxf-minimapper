#!/usr/bin/env Rscript

# combine all coverage.tsv from all samples
# generate an tsv and html table 

# arg[1] is a pattern to list files
# arg[2] is workflow id
library(vroom)
library(dplyr)
library(DT)
library(stringr)
library(sparkline)

arg <- commandArgs(trailingOnly = T)
tsvfiles <- list.files(pattern = arg[1], full.names = T)

cols <- c(
    "qname", "rname", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq", "cov_depth", "cons_qual" 
)
# return a composite sparkline for coverage depth and consensus quality
# x,y are the chr column with values
make_spk <- function(x,y) {
  sp1 <- sparkline(
    str_split(x, "\\|") %>% unlist() %>% head(-1),
    width = 400, height = 40, lineColor = '#2b83ba', fillColor = "#d6eaf8", lineWidth = 2.0, chartRangeMin = 0,
    tooltipFormat = "<span style='color: {{color}}'>&#9679;</span> {{prefix}}pos: {{x}} depth: {{y}} {{suffix}}</span>"
  )
  sp2 <- sparkline(
    str_split(y, "\\|") %>% unlist() %>% head(-1),
    width = 400, height = 40, lineColor = '#d7191c', fillColor = FALSE, lineWidth = 2.0, chartRangeMin = 0,
    tooltipFormat = "<span style='color: {{color}}'>&#9679;</span> {{prefix}}pos: {{x}} cons_q: {{y}} {{suffix}}</span>"
  )
  sl <- spk_composite(sp1, sp2)
  as.character(htmltools::as.tags(sl))
}

df <- vroom(file = tsvfiles, delim = "\t", col_names = cols) #col_select = -c('startpos', 'endpos'))
refsize <- df$endpos[1]

df2 <- df %>% 
  dplyr::select(-c('startpos', 'endpos')) %>%
  group_by(qname) %>%
  mutate(
    cov_depth = make_spk(cov_depth, cons_qual)
    # cov_depth = spk_chr(
    #   str_split(cov_depth, "\\|") %>% unlist() %>% 
    #     #str_replace(pattern = ".*:0$", replacement = ":null") %>% # replace 0 with null so that no line is drawn
    #     head(-1), # remove last element because it is "". This comes from `tr "\n" "|"` 
    #   width = 400, height = 40, lineColor = '#566573', fillColor = '#e5f5e0', lineWidth = 2.0, chartRangeMin = 0,
    #   tooltipFormat = "<span style='color: {{color}}'>&#9679;</span> {{prefix}}pos: {{x}} depth: {{y}} {{suffix}}</span>"
    #   )
    ) %>%
  dplyr::select(-c(cons_qual))
  #dplyr::relocate(coverage, .after = meanmapq)

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
    dplyr::arrange(df2, stringi::stri_rank(qname, opts_collator = locale)),
    class = 'row-border hover',
    # caption = paste0("Run name: ", arg[2], " | Time: ", format.POSIXct(Sys.time())),
    caption = htmltools::tags$caption(
      style = 'caption-side: bottom; text-align: left; color: grey;',
      htmltools::HTML(
        "Ref size: <i>", refsize, 
        "bp </i><br/>Run id:&nbsp&nbsp <i>", arg[2],
        "</i><br/>Date:&nbsp&nbsp&nbsp&nbsp&nbsp <i>", format.POSIXct(Sys.time(), format = "%Y-%m-%d")
        )
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
  DT::formatRound('meandepth', digits = 0, mark = "") %>%
  DT::formatStyle('coverage', color = styleInterval(c(50, 80), c('#e74c3c', '#f5b041', '#1e8449'))) %>%
  DT::formatStyle(c('meanbaseq', 'meandepth'), color = styleInterval(c(20, 25), c('#e74c3c', '#f5b041', '#1e8449'))) %>%
  DT::formatStyle('meanmapq', color = styleInterval(c(40, 50), c('#e74c3c', '#f5b041', '#1e8449'))) %>%
  spk_add_deps()

#write.csv(df, file = '00-alignment-summary.tsv', sep = '\t', row.names = F, col.names = T)
DT::saveWidget(finaltable, file = '00-alignment-summary.html', title = "minimapper-summary")