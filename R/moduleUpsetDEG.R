# modules/upsetDegUI.R
#' Sidebar
#' @noRd
# modules/upsetDegUI.R
# ===================== UI =====================
upsetDegUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Datasets
    selectInput(
      ns("pick_datasets"),
      label   = "Datasets",
      choices = NULL,
      multiple = TRUE,
      selectize = TRUE
    ),

    # Up / Down / Both
    radioButtons(
      ns("direction"), "Direction",
      choices = c("Up", "Down", "Both"),
      inline = TRUE, selected = "Up"
    ),

    # Intersection rule (default to ≥ N since ALL is empty in your data)
    radioButtons(
      ns("mode"), "Select genes by:",
      choices = c(
        "Present in ≥ N datasets" = "atleast",
        "Common across ALL"       = "all",
        "Custom combination"      = "custom"
      ),
      inline = FALSE, selected = "atleast"
    ),

    # Mode-specific controls: rendered by server
    uiOutput(ns("mode_controls")),

    # Add after uiOutput(ns("mode_controls"))
    numericInput(ns("top_k"), "Show top intersections", value = 30, min = 5, max = 200, step = 5),
    numericInput(ns("min_inter_size"), "Hide intersections smaller than", value = 10, min = 1, max = 1000, step = 1),

    actionButton(ns("go"), "Update plot", class = "btn btn-primary"),

    # Helpful summary (max achievable, current N, result size)
    div(style="margin-top:6px; margin-bottom:6px; font-size: 12px; color:#444;",
        textOutput(ns("summary_text"), inline = TRUE)
    ),

    # UpSet plot
    shinycssloaders::withSpinner(
      upsetjs::upsetjsOutput(ns("upset"), height = "680px"),
      type = 4,  color = "#005249"
    ),

    tags$hr(),
    fluidRow(
      column(6, h4("Selected genes")),
      column(6, div(style="text-align:right; margin-top:8px;",
                    downloadButton(ns("download_genes"), "Download")))
    ),
    shinycssloaders::withSpinner(
      DT::dataTableOutput(ns("genes_table")),
      type = 4, color = "#005249"
    )
  )
}

# modules/upsetDegServer.R
#' @noRd
upsetDegServer <- function(id, deg_long) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    `%||%` <- function(a, b) if (is.null(a)) b else a

    cat("[DEBUG] upsetDegServer loaded\n")

    # 1) Populate dataset choices (no default selection)
    observe({
      ds <- sort(unique(deg_long$dataset))
      updateSelectInput(session, "pick_datasets", choices = ds, selected = character(0))
    })

    # 2) Filtered long table (cached)
    filtered_long <- reactive({
      req(length(input$pick_datasets) > 0)
      dl <- deg_long[deg_long$dataset %in% input$pick_datasets, , drop = FALSE]
      dir_sel <- input$direction
      if (!is.null(dir_sel) && dir_sel %in% c("Up","Down")) {
        dl <- dl[dl$direction == dir_sel, , drop = FALSE]
      }
      cat("\n[DEBUG] filtered_long(): datasets =", paste(input$pick_datasets, collapse=","),
          " direction =", dir_sel, " rows =", nrow(dl), "\n")
      dl
    }) |>
      bindCache(input$pick_datasets, input$direction)

    # 3) Build sets list (split first, then unique) — cached
    sets_list_raw <- reactive({
      dl <- filtered_long()
      lapply(split(dl$gene, dl$dataset), unique)
    }) |>
      bindCache(input$pick_datasets, input$direction)

    # ---- FAST: build counts via bitmasks (dimension-safe) -------------------
    build_upset_expression <- function(sets, top_k = 30, min_size = 10) {
      if (!length(sets)) return(list())

      all_genes <- unique(unlist(sets, use.names = FALSE))
      if (!length(all_genes)) return(list())

      gid <- seq_along(all_genes); names(gid) <- all_genes
      ds_names <- names(sets)
      S <- length(ds_names)
      if (S == 0L) return(list())

      # Incidence matrix (genes x datasets)
      M <- matrix(FALSE, nrow = length(all_genes), ncol = S)
      for (j in seq_len(S)) {
        g <- sets[[j]]
        idx <- gid[g]
        idx <- idx[!is.na(idx)]
        if (length(idx)) M[idx, j] <- TRUE
      }
      if (!nrow(M) || !ncol(M)) return(list())

      Mi <- matrix(as.integer(M), nrow = nrow(M), ncol = ncol(M))
      w  <- bitwShiftL(1L, 0:(S-1L))
      mask <- as.vector(Mi %*% matrix(as.integer(w), ncol = 1L))

      tab <- tapply(mask, mask, length)
      if (is.null(tab)) return(list())
      tab <- sort(tab[names(tab) != "0"], decreasing = TRUE)
      if (!length(tab)) return(list())

      tab <- tab[tab >= as.integer(min_size)]
      if (!length(tab)) return(list())
      if (length(tab) > top_k) tab <- tab[seq_len(top_k)]

      key_from_mask <- function(m) {
        on <- which(bitwAnd(as.integer(m), w) != 0L)
        paste(ds_names[on], collapse = "&")
      }
      keys <- vapply(as.integer(names(tab)), key_from_mask, character(1))

      out <- as.list(as.integer(tab))
      names(out) <- keys
      out
    }
    # ------------------------------------------------------------------------

    # Helper: compute max achievable overlap quickly (dimension-safe)
    compute_max_overlap <- function(sets) {
      if (!length(sets)) return(0L)
      all_genes <- unique(unlist(sets, use.names = FALSE))
      if (!length(all_genes)) return(0L)

      gid <- seq_along(all_genes); names(gid) <- all_genes
      S <- length(sets)
      M <- matrix(FALSE, nrow = length(all_genes), ncol = S)
      for (j in seq_len(S)) {
        g <- sets[[j]]
        idx <- gid[g]
        idx <- idx[!is.na(idx)]
        if (length(idx)) M[idx, j] <- TRUE
      }
      if (!nrow(M)) return(0L)
      as.integer(max(rowSums(M)))
    }

    # 4) Event reactives (triggered by Update button)
    sets_list_ev <- eventReactive(input$go, {
      s <- sets_list_raw()
      cat("[DEBUG] sets_list(): #sets =", length(s), "\n")
      s
    }, ignoreInit = TRUE)

    max_overlap_ev <- eventReactive(input$go, {
      s <- sets_list_raw()
      m <- compute_max_overlap(s)
      cat("[DEBUG] max_overlap():", m, "\n")
      m
    }, ignoreInit = TRUE)

    # 5) Mode-specific controls (custom include/exclude)
    output$mode_controls <- renderUI({
      mode <- input$mode %||% "atleast"
      if (identical(mode, "custom")) {
        ds <- names(sets_list_raw())
        fluidRow(
          column(6, selectInput(ns("must_include"), "Must include (present in ALL):",
                                choices = ds, multiple = TRUE, selectize = FALSE, size = 8)),
          column(6, selectInput(ns("must_exclude"), "Must exclude (present in ANY):",
                                choices = ds, multiple = TRUE, selectize = FALSE, size = 8))
        )
      } else NULL
    })

    # Default N (auto)
    get_safe_n <- reactive({
      if (length(sets_list_raw()) <= 1L) 1L else 2L
    })

    # 6) Selected genes (on click)
    selected_genes_ev <- eventReactive(input$go, {
      s <- sets_list_raw()
      if (!length(s)) return(character(0))
      mode <- input$mode %||% "atleast"

      if (mode == "all") {
        if (length(s) == 1) s[[1]] else Reduce(intersect, s)
      } else if (mode == "atleast") {
        n <- get_safe_n()
        all_genes <- unique(unlist(s, use.names = FALSE))
        if (!length(all_genes)) return(character(0))
        sets_env <- lapply(s, function(v) { e <- new.env(hash=TRUE, parent=emptyenv()); for (g in v) assign(g, TRUE, envir=e); e })
        counts <- vapply(all_genes, function(g)
          sum(vapply(sets_env, function(e) exists(g, envir=e, inherits=FALSE), logical(1))),
          integer(1))
        all_genes[counts >= n]
      } else {
        inc <- input$must_include; exc <- input$must_exclude
        base <- if (length(inc)) Reduce(intersect, s[inc]) else unique(unlist(s, use.names = FALSE))
        if (length(exc)) base <- setdiff(base, unique(unlist(s[exc], use.names = FALSE)))
        base
      }
    }, ignoreInit = TRUE)

    # 7) Summary text (on click)
    output$summary_text <- renderText({
      req(input$go)
      s    <- sets_list_raw()
      mode <- input$mode %||% "atleast"
      nset <- length(s); if (!nset) return("No datasets selected.")
      maxAch <- as.integer(max_overlap_ev() %||% 0L)
      sel    <- length(selected_genes_ev())

      if (mode == "all") {
        if (sel == 0) return(paste0("Mode: ALL. Datasets: ", nset,
                                    ". Max achievable overlap: ", maxAch,
                                    ". No genes are shared across all selected datasets."))
        paste0("Mode: ALL. Datasets: ", nset,
               ". Max achievable overlap: ", maxAch,
               ". Genes shared across all: ", sel, ".")
      } else if (mode == "atleast") {
        n <- get_safe_n()
        if (isTRUE(n > maxAch)) {
          paste0("Mode: ≥N (auto N=", n, "). Datasets: ", nset,
                 ". Max achievable overlap = ", maxAch,
                 ". No genes meet ≥", n, " with current selections.")
        } else {
          paste0("Mode: ≥N (auto N=", n, "). Datasets: ", nset,
                 ". Max achievable overlap = ", maxAch,
                 ". Genes meeting ≥", n, ": ", sel, ".")
        }
      } else {
        paste0("Mode: custom. Datasets: ", nset, ". Genes selected: ", sel, ".")
      }
    })

    # 8) UpSet plot (on click)
    output$upset <- upsetjs::renderUpsetjs({
      req(input$go)
      s <- sets_list_ev()
      if (!length(s)) return(upsetjs::upsetjs())

      expr <- build_upset_expression(
        sets      = s,
        top_k     = as.integer(input$top_k %||% 30),
        min_size  = as.integer(input$min_inter_size %||% 10)
      )

      if (!length(expr)) {
        upsetjs::upsetjs() |> upsetjs::fromList(s)
      } else {
        upsetjs::upsetjs() |> upsetjs::fromExpression(expr)
      }
    })

    # 9) Table + download (on click)
    output$genes_table <- DT::renderDataTable({
      req(input$go)
      g <- sort(unique(selected_genes_ev()))
      if (!length(g)) {
        DT::datatable(data.frame(gene = character(0)),
                      options = list(dom = 't',
                                     language = list(emptyTable = "No genes match the current rule.")),
                      rownames = FALSE)
      } else {
        DT::datatable(data.frame(gene = g),
                      options = list(pageLength = 12, scrollX = TRUE),
                      rownames = FALSE)
      }
    })

    output$download_genes <- downloadHandler(
      filename = function() {
        mode <- input$mode %||% "atleast"
        nval <- if (identical(mode, "atleast")) paste0("_N", get_safe_n()) else ""
        paste0("intersection_genes_", (input$direction %||% "Both"), "_", mode, nval, ".txt")
      },
      content = function(file) {
        writeLines(sort(unique(selected_genes_ev())), con = file)
      }
    )
  })
}
