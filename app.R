library(shiny)
library(bslib)
library(ggplot2)

# ── Karyotype definitions ─────────────────────────────────────────────────────
# Each entry: minor allele copies (a) and total copies (t) in the tumor cell
karyotypes <- list(
  "A"    = list(a = 0, t = 1, desc = "Hemizygous deletion"),
  "AA"   = list(a = 0, t = 2, desc = "Copy-neutral LOH"),
  "AB"   = list(a = 1, t = 2, desc = "Diploid (same as normal) — indeterminate"),
  "AAB"  = list(a = 1, t = 3, desc = "Copy Gain"),
  "AAA"  = list(a = 0, t = 3, desc = "Triploid, all major allele"),
  "AABB" = list(a = 2, t = 4, desc = "Tetraploid balanced (same as normal) — indeterminate"),
  "AAAB" = list(a = 1, t = 4, desc = "Tetraploid, 1 minor allele"),
  "AAAA" = list(a = 0, t = 4, desc = "Tetraploid, all major allele")
)

is_indeterminate <- function(k) k$a / k$t == 0.5

# p = (1 - 2m) / (m*(t-2) + (1-a))
calc_purity <- function(maf, a, t) {
  denom <- maf * (t - 2) + (1 - a)
  if (abs(denom) < 1e-10) return(NA_real_)
  p <- (1 - 2 * maf) / denom
  if (p < -1e-6 || p > 1 + 1e-6) return(NA_real_)
  round(min(max(p, 0), 1), 4)
}

# Observed MAF given purity for a karyotype
maf_from_purity <- function(p, a, t) {
  (p * a + (1 - p) * 1) / (p * t + (1 - p) * 2)
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- page_sidebar(
  title = "Tumor Purity Calculator",
  theme = bs_theme(
    bootswatch = "flatly",
    primary    = "#2c7bb6"
  ),

  sidebar = sidebar(
    width = 300,

    h5("Inputs", class = "text-primary fw-semibold"),

    selectInput(
      "karyotype", "Tumor event karyotype",
      choices  = names(karyotypes),
      selected = "AA"
    ),

    uiOutput("karyotype_info"),

    hr(),

    numericInput(
      "maf", "Observed minor allele fraction (MAF)",
      value = 0.3, min = 0, max = 0.5, step = 0.01
    ),

    uiOutput("maf_hint"),

    hr(),

    actionButton(
      "calc", "Calculate purity",
      class = "btn-primary w-100",
      icon  = icon("calculator")
    ),

    hr(),

    p(class = "text-muted small",
      "Assumes the CNV event is clonal (100% clonality). Normal cells are diploid AB.",
      "Formula: p = (1 − 2m) / (m·(t − 2) + (1 − a))")
  ),

  # ── Main panel ──────────────────────────────────────────────────────────────
  layout_columns(
    col_widths = c(4, 8),

    # Result card
    card(
      card_header("Estimated Tumor Purity"),
      card_body(
        uiOutput("result_display")
      )
    ),

    # Plot card
    card(
      card_header("MAF vs. Tumor Purity Curve"),
      card_body(
        plotOutput("purity_curve", height = "340px")
      )
    )
  ),

  # Details card
  card(
    card_header("Calculation Details"),
    card_body(
      tableOutput("detail_table")
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  k <- reactive(karyotypes[[input$karyotype]])

  # Dynamic karyotype description
  output$karyotype_info <- renderUI({
    info <- k()
    tags$p(
      class = "text-muted small mt-1",
      icon("info-circle"),
      info$desc,
      tags$br(),
      sprintf("Copies: %d minor + %d major = %d total",
              info$a, info$t - info$a, info$t)
    )
  })

  # MAF hint: valid range
  output$maf_hint <- renderUI({
    ki <- k()
    if (is_indeterminate(ki)) {
      tags$p(class = "text-warning small",
             icon("exclamation-triangle"),
             "This karyotype cannot be distinguished from normal. Purity is indeterminate.")
    } else {
      lo <- ki$a / ki$t
      tags$p(class = "text-muted small",
             sprintf("Valid range for this karyotype: [%.3f, 0.500]", lo))
    }
  })

  # Computed purity (reactive, triggered by button)
  purity <- eventReactive(input$calc, {
    ki  <- k()
    maf <- input$maf

    if (is_indeterminate(ki)) return(list(p = NA, error = "Indeterminate karyotype"))
    if (is.na(maf) || maf < 0 || maf > 0.5) return(list(p = NA, error = "MAF must be between 0 and 0.5"))

    lo <- ki$a / ki$t
    if (maf < lo - 1e-6) {
      return(list(p = NA,
                  error = sprintf("MAF %.4f is below the minimum possible (%.4f) for this karyotype.", maf, lo)))
    }

    p <- calc_purity(maf, ki$a, ki$t)
    if (is.na(p)) return(list(p = NA, error = "No valid purity solution for these inputs."))
    list(p = p, error = NULL)
  }, ignoreNULL = FALSE)

  # Result display
  output$result_display <- renderUI({
    res <- purity()
    if (!is.null(res$error)) {
      tagList(
        tags$div(class = "alert alert-warning", icon("exclamation-triangle"), " ", res$error)
      )
    } else {
      pct <- res$p * 100
      tagList(
        tags$div(
          class = "text-center py-3",
          tags$span(
            class = "display-4 fw-bold text-primary",
            sprintf("%.1f%%", pct)
          ),
          tags$p(class = "text-muted mt-2",
                 sprintf("(p = %.4f)", res$p))
        )
      )
    }
  })

  # Detail table
  output$detail_table <- renderTable({
    ki  <- k()
    maf <- input$maf
    res <- purity()

    df <- data.frame(
      Parameter = c(
        "Karyotype",
        "Minor copies in tumor (a)",
        "Total copies in tumor (t)",
        "Tumor MAF (a/t)",
        "Normal MAF",
        "Observed MAF (input)",
        "Valid MAF range",
        "Tumor purity (p)"
      ),
      Value = c(
        input$karyotype,
        ki$a,
        ki$t,
        sprintf("%.4f", ki$a / ki$t),
        "0.5000",
        sprintf("%.4f", maf),
        if (is_indeterminate(ki)) "N/A"
        else sprintf("[%.4f, 0.5000]", ki$a / ki$t),
        if (!is.null(res$error)) paste("Error:", res$error)
        else sprintf("%.4f  (%.1f%%)", res$p, res$p * 100)
      ),
      stringsAsFactors = FALSE
    )
    df
  }, striped = TRUE, hover = TRUE, bordered = TRUE, width = "100%")

  # MAF vs purity curve
  output$purity_curve <- renderPlot({
    ki  <- k()
    maf <- input$maf
    res <- purity()

    if (is_indeterminate(ki)) {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, size = 5, color = "gray50",
                 label = "This karyotype is indeterminate\n(tumor MAF = normal MAF = 0.5)") +
        theme_void()
    } else {
      purities <- seq(0.001, 1, length.out = 500)
      mafs     <- maf_from_purity(purities, ki$a, ki$t)
      df       <- data.frame(purity = purities, maf = mafs)

      gg <- ggplot(df, aes(x = purity, y = maf)) +
        geom_line(color = "#2c7bb6", linewidth = 1.3) +
        scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                           limits = c(0, 1), expand = expansion(mult = 0.01)) +
        scale_y_continuous(limits = c(0, 0.5), expand = expansion(mult = 0.02)) +
        labs(
          x = "Tumor purity",
          y = "Observed MAF",
          title = sprintf("Karyotype: %s  |  a = %d, t = %d",
                          input$karyotype, ki$a, ki$t)
        ) +
        theme_minimal(base_size = 13) +
        theme(plot.title = element_text(face = "bold", size = 13))

      # Mark the current input MAF as a horizontal line
      gg <- gg +
        geom_hline(yintercept = maf, linetype = "dashed",
                   color = "gray60", linewidth = 0.8)

      # Mark the solution point if valid
      if (!is.null(res) && is.null(res$error)) {
        gg <- gg +
          geom_point(data = data.frame(purity = res$p, maf = maf),
                     aes(x = purity, y = maf),
                     color = "#d7191c", size = 4, shape = 16) +
          geom_vline(xintercept = res$p, linetype = "dashed",
                     color = "#d7191c", linewidth = 0.8) +
          annotate("text",
                   x     = min(res$p + 0.03, 0.97),
                   y     = maf + 0.025,
                   label = sprintf("p = %.1f%%", res$p * 100),
                   color = "#d7191c", fontface = "bold", size = 4.5,
                   hjust = 0)
      }
      gg
    }
  }, res = 110)
}

shinyApp(ui, server)
