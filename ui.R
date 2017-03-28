library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)





source('PowerCalc_Rare.r')


# Define UI for slider demo application
tagList(
  navbarPage(
    theme = shinythemes::shinytheme("cerulean"),
    "AATAP: Aggregated Association test Approximate Power",
    tabPanel("Overview",
             fluidPage(
               # sidebarPanel(
               #   p("AATAP is designed by Andriy Derkach, Haoyu Zhang and Nilanjan Chatterjee", style = "font-family: 'times'; font-si16pt")
               # ),
               mainPanel(
                 h3("Introduction"),
                 p("These are simplified power calculations for commonly used rare variant tests such as SKAT, C-alpha, Hotelling, Burden [Wu et al. 2011; Neale et al. 2011; Derkach et al. 2014].  Application estimates average power of at test in genome-wide study by using minimum number of parameters:
", style = "font-family: 'times'; font-si20pt"),
                 h3("Essential Parameters"),
                 p("1)	EV=% of variation in a trait explained by loci (gene); 

2) alpha = level of the test; 

3) Number of cases/Number of controls in case-control study; 

4) Total Sample Size in continuous trait study. 
", style = "font-family: 'times'; font-si20pt"),
                 h3("Optional Parameters
"),
                 p("1)	Range of EV= range of percentages of variation in a trait explained by loci; 2) Proportion of Causal = proportion of causal variants in a locus; 3) Number of Variants = number of causal variants in a locus.
", style = "font-family: 'times'; font-si20pt"),
                 h3("Assumptions"),
                 p("Currently power calculations are based on public data for 60,000 individuals from Exome Aggregation Consortium [Lek et al. 2016]. It may underestimate number of variants per gene in whole-genome study (it overestimates an average power). Power calculations estimate an average power under three relationships: 1) there is no relationship between MAF and % of variations explained by a variant; 2) there is no relationship between MAF and effect size (log-OR) and 3) effect size (log-OR) is proportional to log10(MAF). 
", style = "font-family: 'times'; font-si20pt"),
                 h3("Contributor"),
                 p("AATAP is designed by Andriy Derkach, Haoyu Zhang and Nilanjan Chatterjee
", style = "font-family: 'times'; font-si16pt")
                 ,width = 15)
             )),
    tabPanel("Case Control",
            sidebarPanel( 
              selectInput(
                "method", "Method",
                c("SKAT","Calpha","Hotelling","Burden Test")
                ),
             
              numericInput("ncases", "Number of Cases:", 
                                       value=5000),
                          
                          
                          numericInput("ncont", "Number of Controls:", 
                                       value = 5000),
                          
                          numericInput("Alpha",HTML("&alpha;:"), value=0.0001),
                          
              
                          selectInput("evoption", 
                                      "Do you want to specify a range of EV?",
                                      c("No","Yes")),
              conditionalPanel(
                condition = "input.evoption == 'No'",
                numericInput("EV", "EV: Variance Explained(Percent)",
                             value = 1)),
              conditionalPanel(
                condition = "input.evoption == 'Yes'",
                sliderInput("EV_range","Variance Explained Range(Percent)",
                            min=0.01,max=1,value = c(0.05,0.60),step=0.01)),
              
              # Only show this panel if the plot type is a histogram
             
               conditionalPanel(
               condition="input.method!='Burden Test'",
               numericInput("PC","Proportion of Causal (Optional)",value=NULL)
               ),   
              conditionalPanel(
                condition="input.method=='Burden Test'",
                numericInput("PC_new","Proportion of Causal",value=0.05),
                numericInput("PRC","Proportion of Protective",value=0)
              ),
                         
      # textInput("nameEsseble","Specify One Gene(Optional)"),
       numericInput("JJ","Number of Variants (Optional)",value=NULL),
              
          
      actionButton("update","Update")
            ),
             # Show a table summarizing the values entered
             mainPanel(
               tabsetPanel(
                 tabPanel("Power Density",
                          h4("The result in this page will show up if you choose NOT to specify EV range"),
                          p(HTML("Genetic Architecture I: &beta; independent of MAF"),style = "font-family: 'times'; font-si20pt"),
                          p("Geneic Architecture II : EV independent of MAF",style = "font-family: 'times'; font-si20pt"), 
                          p(HTML("Genetic Architecture III: &beta;~log_10(MAF);"),style = "font-family: 'times'; font-si20pt"),
                          dataTableOutput("values"),
                           plotOutput("plot")
                           ),
               tabPanel("Mean Power Distribution based on different EV",
                        h4("The result in this page will show up if you choose to specify a range of EV"),
                           p(HTML("Genetic Architecture I: &beta; independent of MAF"),style = "font-family: 'times'; font-si20pt"),
                           p("Geneic Architecture II : EV independent of MAF",style = "font-family: 'times'; font-si20pt"), 
                           p(HTML("Genetic Architecture III: &beta;~log_10(MAF);"),style = "font-family: 'times'; font-si20pt"),
                       
                          # dataTableOutput("values"),
                        #dataTableOutput("values_r"),
                        dataTableOutput("values_r"),
                          plotOutput("plot_r")
             )
    )
    )
    ),
    tabPanel("QT",
             sidebarPanel( 
               selectInput(
                 "method2", "Method",
                 c("SKAT","Calpha","Hotelling","Burden Test")
               ),
              
               numericInput("total2", "Total Sample Size:", 
                            value=10000),
               
               numericInput("Alpha2",HTML("&alpha;:"), value=0.0001
               ),
               
               selectInput(
                 "evoption2", "Do you want to specify a range of EV?",
                 c("No","Yes")),
               conditionalPanel(
                 condition = "input.evoption2 == 'No'",
                 numericInput("EV2", "EV: Variance Explained(Percent)",
                              value = 0.5)),
               # Only show this panel if the plot type is a histogram
               conditionalPanel(
                 condition = "input.evoption2 == 'Yes'",
                 sliderInput("EV_range2","Variance Explained Range(Percent)",
                             min=0.01,max=1,value = c(0.05,0.60),step=0.01)),
               conditionalPanel(
                 condition="input.method2!='Burden Test'",
                 numericInput("PC2","Proportion of Causal (Optional)",value=NULL)
                ),
               conditionalPanel(
                 condition="input.method2=='Burden Test'",
                 numericInput("PC_new2","Proportion of Causal",value=0.05),
                 numericInput("PRC2","Proportion of Protective",value=0)
               ),
             #  textInput("nameEsseble2","Specify One Gene(Optional)"),
               numericInput("JJ2","Number of Variants (Optional)",value=NULL),
              
               
               actionButton("update2","Update")
             ),
             mainPanel(
             
               tabsetPanel(
                 tabPanel("Power Density",
                          h4("The result in this page will show up if you choose NOT to specify EV range"),     
                          p(HTML("Genetic Architecture I: &beta; independent of MAF"),style = "font-family: 'times'; font-si20pt"),
                          p("Geneic Architecture II : EV independent of MAF",style = "font-family: 'times'; font-si20pt"), 
                          p(HTML("Genetic Architecture III: &beta;~log_10(MAF);"),style = "font-family: 'times'; font-si20pt"),
                          dataTableOutput("values2"),
                          plotOutput("plot2")
                 ),
                 tabPanel("Mean Power Distribution based on different EV",
                          h4("The result in this page will show up if you choose to specify a range of EV"),
                          p(HTML("Genetic Architecture I: &beta; independent of MAF"),style = "font-family: 'times'; font-si20pt"),
                          p("Geneic Architecture II : EV independent of MAF",style = "font-family: 'times'; font-si20pt"), 
                          p(HTML("Genetic Architecture III: &beta;~log_10(MAF);"),style = "font-family: 'times'; font-si20pt"),
                          dataTableOutput("values2_r"),
                          plotOutput("plot2_r")
                 )
               )
             )
            
    ),
    tabPanel("Citation")
)
  )



