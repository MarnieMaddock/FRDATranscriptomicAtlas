# modules/upsetDegUI.R
#' Sidebar
#' @noRd
upsetDegUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Choose datasets
    selectizeInput(
      ns("pick_datasets"),
      label = "Datasets",
      choices = NULL,
      multiple = TRUE,
      options = list(plugins = list("remove_button"), placeholder = "Select one or more datasets")
    ),

    # Up/Down/Both
    radioButtons(
      ns("direction"), "Direction",
      choices = c("Up", "Down", "Both"),
      inline = TRUE, selected = "Up"
    ),

    # How to select intersection genes
    radioButtons(
      ns("mode"), "Select genes by:",
      choices = c(
        "Common across ALL" = "all",
        "Present in ≥ N datasets" = "atleast",
        "Custom combination" = "custom"
      ),
      inline = FALSE, selected = "all"
    ),

    # Controls that depend on mode (rendered in server)
    uiOutput(ns("mode_controls")),

    # UpSet plot
    upsetjs::upsetjsOutput(ns("upset"), height = "460px"),

    tags$hr(),
    fluidRow(
      column(6, h4("Selected genes")),
      column(6, div(style="text-align:right; margin-top:8px;",
                    downloadButton(ns("download_genes"), "Download")))
    ),
    DT::dataTableOutput(ns("genes_table"))
  )
}

#' Server
#' @noRd
upsetDegServer <- function(id, deg_long) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # 1) Populate dataset choices
    observe({
      ds <- sort(unique(deg_long$dataset))
      updateSelectizeInput(session, "pick_datasets",
                           choices = ds, selected = ds, server = TRUE)
    })

    # 2) Build per-dataset gene sets (after Direction filter)
    sets_list <- reactive({
      req(input$pick_datasets)
      dl <- deg_long[deg_long$dataset %in% input$pick_datasets, , drop = FALSE]

      dir_sel <- input$direction
      if (!is.null(dir_sel) && dir_sel %in% c("Up","Down")) {
        dl <- dl[dl$direction == dir_sel, , drop = FALSE]
      }
      # Both = no extra filter

      # split to a named list: dataset -> unique gene vector
      l <- split(dl$gene, dl$dataset)
      l <- lapply(l, function(v) sort(unique(v)))
      l
    })

    # 3) UpSet plot (exported API only)
    output$upset <- upsetjs::renderUpsetjs({
      s <- sets_list()
      # Render even with 1 set (shows composition); UpSet shines with ≥2
      upsetjs::upsetjs() |> upsetjs::fromList(s)
    })

    # 4) Mode-specific controls (slider or custom include/exclude)
    output$mode_controls <- renderUI({
      ds <- names(sets_list())
      mode <- input$mode %||% "all"
      if (mode == "atleast") {
        max_n <- max(2L, length(ds))
        sliderInput(ns("min_n"), "N (at least)", min = 2, max = max_n, value = min(2, max_n), step = 1)
      } else if (mode == "custom") {
        fluidRow(
          column(6,
                 selectizeInput(ns("must_include"), "Must include (present in ALL):",
                                choices = ds, multiple = TRUE,
                                options = list(plugins = list("remove_button")))),
          column(6,
                 selectizeInput(ns("must_exclude"), "Must exclude (present in ANY):",
                                choices = ds, multiple = TRUE,
                                options = list(plugins = list("remove_button"))))
        )
      } else {
        # mode "all" → no extra controls
        NULL
      }
    })

    # 5) Compute selected genes per mode
    selected_genes <- reactive({
      s <- sets_list()
      if (length(s) == 0) return(character(0))

      mode <- input$mode %||% "all"

      if (mode == "all") {
        # Intersection across ALL chosen datasets
        if (length(s) == 1) return(s[[1]])
        return(Reduce(intersect, s))
      }

      if (mode == "atleast") {
        n <- input$min_n
        if (is.null(n) || is.na(n)) n <- 2L
        n <- max(1L, min(as.integer(n), length(s)))
        all_genes <- unique(unlist(s, use.names = FALSE))
        # Count in how many sets each gene appears
        counts <- vapply(all_genes, function(g) {
          sum(vapply(s, function(vec) g %in% vec, logical(1)))
        }, integer(1))
        return(all_genes[counts >= n])
      }

      # mode == "custom"
      inc <- input$must_include
      exc <- input$must_exclude

      # Start from intersection of 'inc' if provided; else union of all
      base <- if (length(inc)) {
        Reduce(intersect, s[inc])
      } else {
        unique(unlist(s, use.names = FALSE))
      }

      # Remove any gene present in ANY excluded set
      if (length(exc)) {
        ex_union <- unique(unlist(s[exc], use.names = FALSE))
        base <- setdiff(base, ex_union)
      }
      base
    })

    # 6) Show table + download
    output$genes_table <- DT::renderDataTable({
      g <- sort(unique(selected_genes()))
      DT::datatable(
        data.frame(gene = g),
        options = list(pageLength = 12, scrollX = TRUE),
        rownames = FALSE
      )
    })

    output$download_genes <- downloadHandler(
      filename = function() {
        paste0("intersection_genes_", (input$direction %||% "Both"), "_", (input$mode %||% "all"), ".txt")
      },
      content = function(file) {
        writeLines(sort(unique(selected_genes())), con = file)
      }
    )
  })
}

`%||%` <- function(a, b) if (is.null(a)) b else a
