team_tab <- function() {
  tabPanel(
    "Team",
    value = "team",
    div(
      class = "page-container",   # <-- This aligns with navbar & homepage
      
      br(),
      h2("DiCE Development Team"),
      br(),
      
      # --- Sandali ---
      h4("Sandali Lokuge"),
      p(tags$b("Ph.D. Student")),
      p("Department of BioHealth Informatics, Luddy School of Informatics and Computing,"),
      p("Indiana University at Indianapolis, Indianapolis, IN 46202, United States"),
      p(tags$b("Email:"), " sdlokuge@iu.edu"),
      br(),
      
      # --- Sheng Liu ---
      h4("Dr. Sheng Liu"),
      p(tags$b("Research Assistant Professor")),
      p("Department of Medical and Molecular Genetics,"),
      p("Indiana University School of Medicine, Indianapolis, IN 46202, United States"),
      p(tags$b("Email:"), " sliu19@iu.edu"),
      br(),

      
      # --- Jun Wan ---
      h4("Dr. Jun Wan"),
      p(tags$b("Principal Investigator")),
      p("Department of Medical and Molecular Genetics,"),
      p("Indiana University School of Medicine, Indianapolis, IN 46202, United States"),
      p(tags$b("Email:"), " junwan@iu.edu"),
      br()
    )
  )
}
