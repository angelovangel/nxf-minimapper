#!/usr/bin/env Rscript

# combine all coverage.tsv from all samples
# generate an tsv and html table 

# arg[1] is a pattern to list files
# arg[2] is workflow id
# arg[3] is git commit

library(vroom)
library(dplyr)
library(DT)
library(stringr)
library(sparkline)

arg <- commandArgs(trailingOnly = T)
tsvfiles <- list.files(pattern = arg[1], full.names = T)

# 14 columns
cols <- c(
    "qname", "totalreads", "mappedreads", 
    "rname", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq", 
    "cov_depth", "cons_qual" 
)
# return a composite sparkline for coverage depth and consensus quality
# x,y are the chr column with values
make_spk <- function(x,y) {
  sp1 <- sparkline(
    str_split(x, "\\|") %>% unlist() %>% head(-1),
    width = 450, height = 43, lineColor = '#2972b6', fillColor = "#d6eaf8", lineWidth = 2.0, chartRangeMin = 0,
    tooltipFormat = "&#9679; {{prefix}}pos: {{x}} <br><span style='color: {{color}}'>&#9679;</span> depth: {{y}} {{suffix}}</span>"
  )
  sp2 <- sparkline(
    str_split(y, "\\|") %>% unlist() %>% head(-1),
    width = 450, height = 43, lineColor = '#be1a25', fillColor = FALSE, lineWidth = 2.0, chartRangeMin = 0,
    tooltipFormat = "<span style='color: {{color}}'>&#9679;</span> cons_q: {{y}} {{suffix}}</span>"
  )
  sl <- spk_composite(sp1, sp2)
  as.character(htmltools::as.tags(sl))
}

df <- vroom(file = tsvfiles, delim = "\t", col_names = cols) #col_select = -c('startpos', 'endpos'))
refsize <- df$endpos[1]
refname <- df$rname[1]

df2 <- df %>% 
  dplyr::select(-c('startpos', 'endpos', 'rname', 'numreads', 'covbases')) %>%
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
    class = 'row-border',
    colnames = c(
      'Sample',	'Total reads',	'Mapped reads',
      #'Alignments', 
      #'Covered bases', 
      'Coverage%',	'Mean depth',	'Mean base Q', 'Mean map Q', 'Coverage depth'
    ),
    caption = htmltools::tags$caption(
      style = 'caption-side: top; text-align: left; color: grey;',
      htmltools::HTML(
        format.POSIXct(Sys.time(), format = "%Y-%m-%d %H:%M:%S"), 
        "<br/>
        Ref name: <b>", refname, "</b><br/>
        Ref size: &nbsp<b>", refsize, "bp </b><br/>
        Run id:&nbsp&nbsp&nbsp&nbsp<b>", arg[2], "</b><br/><hr/>" 
        #"<hr />"
        )
    ),
    # style = 'bootstrap',
    escape = F, #filter = 'top',
    extensions = 'Buttons', rownames = FALSE,
    options = list(
      searchHighlight = TRUE,
      rowCallback = JS(rowCallback),
      autoWidth = TRUE, pageLength = 125,
      dom = 'Btp',
      paging = FALSE,
      buttons = c('copy', 'csv', 'excel'),
      columnDefs = list(list(width = '300px', targets = 7))
    )
  ) %>% 
  DT::formatRound('coverage', digits = 0) %>%
  DT::formatRound('meandepth', digits = 0, mark = "") %>%
  DT::formatRound('meanbaseq', digits = 0, mark = "") %>%
  DT::formatRound('meanmapq', digits = 0, mark = "") %>%
  DT::formatStyle('coverage', color = styleInterval(c(50, 80), c('#e74c3c', '#f5b041', '#2972b6'))) %>%
  DT::formatStyle(c('meanbaseq', 'meandepth'), color = styleInterval(c(20, 25), c('#e74c3c', '#f5b041', '#2972b6'))) %>%
  DT::formatStyle('meanmapq', color = styleInterval(c(40, 50), c('#e74c3c', '#f5b041', '#2972b6'))) %>%
  spk_add_deps()

# Add Google Fonts link
google_fonts <- htmltools::tags$link(
  href = "https://fonts.googleapis.com/css2?family=Roboto:wght@400;700&display=swap",
  rel = "stylesheet"
)

# Add a modern title with Google Fonts and Nerd Font
intro_text <- htmltools::tags$div(
  style = 'text-align: center; margin-bottom: 20px;',
  htmltools::tags$h1(
    style = 'font-size: 24px; color: #2972b6; margin: 0; font-family: "Roboto", "Hack Nerd Font", monospace;',
    "Alignment Coverage and Quality Metrics" # Example Nerd Font icons
  ),
  htmltools::tags$p(
    style = 'font-size: 15px; color: #7f8c8d; margin-top: 5px; font-family: "Roboto", sans-serif;',
    "Report generated with the ",
    htmltools::tags$a(
      href = "https://github.com/angelovangel/nxf-minimapper",
      target = "_blank",  # Opens the link in a new tab
      style = 'color: #2980b9; text-decoration: none;',
      "nxf-minimapper Nextflow pipeline", htmltools::tags$a("commit:", arg[3])
    )
  )
)

# Prepend the Google Fonts link and title to the table
final_output <- htmlwidgets::prependContent(finaltable, google_fonts, intro_text)

# Save the widget with the combined content
DT::saveWidget(final_output, file = '00-alignment-summary.html', title = "Alignment Coverage Summary", selfcontained = TRUE)
