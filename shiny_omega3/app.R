# ------------------------------------------------------------------------------
# Title: Shiny app script for interactive n-3 FA and IA models
# Author: Ryan Gan
# Date Created: 4/11/2017
# R Version Created Under: 3.3.3
# ------------------------------------------------------------------------------

# load library
library(shiny)
library(tidyverse)
library(broom)

# reading baseline data
baseline_df <- read_csv("baseline.csv") %>% 
  rename(visit_num = visit_) %>% 
  # limit only to variables used in analysis
  select(ia_base, sd_ALA_update, sd_EPA_update, sd_dpa_update, sd_DHA_update,
         sd_epa_dha_update, sd_total_n3_update, age_bi, sex, race_bi, educ, inc,
         eversmoke, cursmoke, py10, se_bi, omega3_bi, multivits_bi, rfneph, crp, 
         ccp2)
 

tab2_matrix <- matrix(nrow = 6, ncol = 5 )
# column names of tables
colnames(tab2_matrix) <- c("n3_FA_RBC", "OR", "Lower_95_CI",
                           "Upper_95_CI", "p_value") 
# variables to loop through
n3_vars <- c("sd_ALA_update", "sd_EPA_update", "sd_dpa_update", 
             "sd_DHA_update", "sd_epa_dha_update", "sd_total_n3_update")
# fill first column of matrix with n3 vars
tab2_matrix[,1] <- c("ALA", "EPA", "DPA", "DHA", "EPA+DHA",
                     "Total n-3 FA")


# Define server logic required to draw a histogram
server <- function(input, output) {

  or_data_input <- reactive({
      # Iteration of logistic regression through n-3 FAs
      data.table(for(i in 1:length(n3_vars)){
        # general logistc regression; using tidy from broom to output 
        # coeffecients
        log_reg_mod <- tidy(glm(as.formula(paste0("ia_base ~", 
                          paste(n3_vars[i], "omega3_bi", "py10", "rfneph", 
                                "crp", "se_bi", sep = "+"))),
                          data = baseline_df, family = "binomial"(link="logit"))) 
        
        # calculate odds ratio and 95% CI
        odds_ratio <- log_reg_mod %>% filter(term == n3_vars[i]) %>% 
          # estimate lower and upper 95% CI bounds
          mutate(lower_bound = estimate - (1.96*std.error),
                 upper_bound = estimate + (1.96*std.error)) %>% 
          select(term, estimate, lower_bound, upper_bound, p.value)
      
        # fill ith row of table 2
        tab2_matrix[i,2:4] <- round(as.numeric(exp(odds_ratio[1, 2:4])),2)
        tab2_matrix[i,5] <- round(as.numeric(odds_ratio[1,5]),2)
        
      }
      )
      }) # end reactive

  
   output$plot <- renderPlot({
     # ggplot of odds ratios
      ggplot(or_data_input(), aes(x = n3_FA_RBC, y = OR)) +
        geom_point(size = 2) +
        geom_errorbar(aes(ymin = Lower_95_CI, ymax = Upper_95_CI),
                      width = 0.15) +
        scale_y_log10(breaks = c(0.01, 0.1, 0.25, 0.5, 0.75, 1, 2.5, 5)) +
        geom_hline(yintercept = 1, linetype=2) +
        ylab('Odds Ratio for 1 Standard Deviation Increase') +
        xlab('Omega-3 Fatty Acid') +
        # plot theme
        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
   })
}


# Define UI for application that draws a histogram
ui <- fluidPage(

   # Application title
   titlePanel("Association Omega-3 Biomarkers and Inflammatory Arthritis"),
   
   plotOutput("plot")
   
)

# Run the application 
shinyApp(ui = ui, server = server)

