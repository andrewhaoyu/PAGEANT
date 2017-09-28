library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)

grid_transform = function(grid){
  if(grid=="quick"){
    return(20)
  }else if(grid=="adequate"){
    return(50)
  }else{
    return(100)
  }
}

busyIndicator <- function (text = "Loading..", img = "busyIndicator/ajaxloaderq.gif", wait = 3) {
  tagList(singleton(tags$head(tags$link(rel = "stylesheet", 
                                        type = "text/css", href = "busyIndicator/busyIndicator.css"))), 
          div(class = "shinysky-busy-indicator", p(text), img(src = img)), 
          tags$script(sprintf("\tsetInterval(function(){\n  \t\t \t if ($('html').hasClass('shiny-busy')) {\n  \t\t    setTimeout(function() {\n  \t\t      if ($('html').hasClass('shiny-busy')) {\n  \t\t        $('div.shinysky-busy-indicator').show()\n  \t\t      }\n  \t\t    }, %d)  \t\t    \n  \t\t  } else {\n  \t\t    $('div.shinysky-busy-indicator').hide()\n  \t\t  }\n  \t\t},100)\n  \t\t", 
                              wait)))
}
source('PowerCalc_Rare.r')
source('calc_whole.r')


# Define UI for slider demo application
tagList(
  navbarPage(
    theme = shinythemes::shinytheme("cerulean"),
    "PAGEANT",
    tabPanel("Overview",
             fluidPage(
               # sidebarPanel(
               #   p("AATAP is designed by Andriy Derkach, Haoyu Zhang and Nilanjan Chatterjee", style = "font-family: 'times'; font-si16pt")
               # ),
               mainPanel(
                 h3("PAGEANT: Power Analysis for GEnetic AssociatioN Tests"),
                 h3("Introduction"),
                 p("The application allows rapid power analysis for a variety of genetic association tests by specification of a few key parameters [Derkach et al. 2017]. Power calculations can be done at the level of a single variant for simple trend test, at the level of a genes/regions for various complex aggregated tests [Neale et al. 2011, Derkach et al. 2013, Wu et al. 2011, Madsen and Browning 2009] and at the level of the whole genome for the assessment of overall yield of a study. The calculations currently uses underlying distribution of gene sizes and minor allele frequencies of variants observed in the in the public data for 60,000 individuals from Exome Aggregation Consortium [Lek et al. 2016] 
", style = "font-family: 'times'; font-si20pt"),
                 h3("Power for association test at the level of a single variant or a single gene/region
"),
                 h4("Essential Input Parameters "),
                 p("1)	EV: For continuous trait, EV represents % of phenotypic variance explained the variants in a gene (or by a single variant for single-variant tests). For binary trait, EV represents sum of squares of log-odds-ratios (in standardized unit) associated with the variants in a gene (or that for a single variant). For example, if power calculation is desired for a locus which may have 5 causal variants each with an OR=1.2, then EV should be set to \\(5× (ln(1.2)^2)=0.17\\);", style = "font-family: 'times'; font-si20pt"),
                 p("2) \\(\\alpha\\) = level of the test;", style = "font-family: 'times'; font-si20pt"),
                 p("3) Sample size: Total sample size for a continuous trait or the number of cases and number of controls for a case-control study;", style = "font-family: 'times'; font-si20pt"),
                 h4("Optional Parameters
"),
                 p("1)	Total number of variants (\\(J\\)): The total number of variants under study within a gene/region. This is a key parameter in power calculation for the gene-level tests and when it’s not specified the application evaluates distribution of power according to distribution observed for in the ExAC database.With J=1, gene based power calculation simplifies that for a single variant.", style = "font-family: 'times'; font-si20pt"),
                 withMathJax(),
                 p("2)	Proportion of causal variants (\\(J_c/J\\)) : Assumed proportion of causal variants in a locus as a ratio to the total number of variants. This parameter is required for burden test and a more accurate second-order approximation of the variance component test. For burden tests, it’s assumed that all causal variants are either deleterious or protective and by default proportion of causal variants is set to 0.2.", style = "font-family: 'times'; font-si20pt"),
                 p("3)Range of EV: Instead of a single EV, the user can specify a range of EV over which power calculation is desired 
                   ", style = "font-family: 'times'; font-si20pt"),
                 h4("Output"),
p("The application conducts power analysis under three different models for genetic architecture assuming (I) MAF is independent of EV (II) MAF is independent of genetic effects measured in the unit of per copy of an allele (\\(\\beta^2\\)=EV/(2MAF(1-MAF)) and (III) MAF is negatively correlated with genetic effect through the function \\(\\beta=-\\log_{10}\\)(MAF).  When a single EV is specified, for each genetic architecture, it returns a distribution of power and key summary measures (mean, median,  25th and 75th percentiles). This distribution corresponds to uncertainty association with various additional parameters, such as the number of variants within a gene and minor allele frequencies. If a range of EV is specified, plots and table for average power over the range of specified EV is returned.

", style = "font-family: 'times'; font-si20pt"),
                 h3("Genome-level power calculation"),
                 h4("Essential Input parameter:"),
                 p("1) M: Hypothesized number of underlying causal loci (or variants if analysis to be done based on single variant test)", style = "font-family: 'times'; font-si20pt"),
                 p("2) GEV: Total EV explained by  loci in genome -wide study (see definition of EV above)", style = "font-family: 'times'; font-si20pt"),
                 

                 h4("Optional Input Parameter"),
                 p(" 1) m: The number of causal loci for which probability of discovery to be calculated (see output).For example if it is assumed M=100, then a user may want to find out what is the probability of discovery of m=5 or less number of loci from a given study of a specific sample size.", style = "font-family: 'times'; font-si20pt"),
p("2) Computational Complexity: The number of  models and iterations used to estimate range of expected number of discoveries and probabilities (see output). There are three options: Quick, Adequate and Complete.  Quick option runs genome-wide calculations  within 3 minutes and  provides rough estimates. Adequate option runs genome-wide calculations  within 6 minutes and provides more accurate estimates. Lastly, Complete option runs genome-wide calculations  within 15 minutes and  provides very accurate estimates. ", style = "font-family: 'times'; font-si20pt"),
                 h4("Output"),
                 p("Expected number of discoveries: The application returns expected number of discoveries where the expectation is calculated across the M loci accounting for uncertainty associated with distribution of number of variants per locus (J), allele frequencies and the distributions of EVs the loci explains. Currently, the distribution of J and MAF in these calculations are obtained from those observed in the ExAC database. In addition, it is assumed the effect size distribution follows a L-shaped gamma distribution with mean specified as \\(\\mu=\\)GEV/M.  The application calculates a range of expected number of discoveries based on the range of the dispersion parameter of the underlying gamma distribution for the effect size distribution and the corresponding maximum and minimum values are returned.  
           Probability of discoveries: This returns maximum and minimum probability of a certain number of discoveries (m) for values of specified by the user. ", style = "font-family: 'times'; font-si20pt"),
                 h3("Additional notes"),
                 p("Currently power calculations are based on distribution of number of variants per gene and minor allele frequencies observed in the public data for 60,000 individuals from Exome Aggregation Consortium [Lek et al. 2016]. It will underestimate the total number of variants per gene/region for whole-genome study and may overestimate power for gene-level tests. ", style = "font-family: 'times'; font-si16pt"),
                 h3("Contributor"),
                 p("PAGEANT is designed by Andriy Derkach, Haoyu Zhang and Nilanjan Chatterjee
", style = "font-family: 'times'; font-si16pt")
                 ,width = 15)
             )),
    tabPanel("Case Control",
             sidebarPanel(
               selectInput("TypeofCalculation", 
                           "Type of calculation",
                           c("Power Calculation","Sample Size Calculation")),
               
               conditionalPanel(
                 condition = "input.TypeofCalculation=='Power Calculation'",
                 selectInput("SNPoption", 
                             "Type of power calculation",
                             c("Single Gene","Single SNP","Whole Genome")),
                 
                 numericInput("ncases", "Number of Cases:", 
                              value=5000),
                 
                 numericInput("ncont", "Number of Controls:", 
                              value = 5000),
                 
                 
                 
                 conditionalPanel(
                   condition = "input.SNPoption=='Whole Genome'",
                   numericInput("Alpha_whole",HTML("&alpha;:"), value=2.5e-06),
                   numericInput("GEV", "GEV: Genome-wide Variance Explained(Percent)",
                                value = 20),
                   selectInput("grid",
                               " Computational Complexity (optional)",
                               c("Quick","Adequate","Complete")),
                   selectInput(
                     "method_whole", "Method",
                     c("SKAT","Calpha","Hotelling","Burden Test")
                   ),
                   numericInput("K","Number of causal loci",value=1000),
                   numericInput("m","Number of discoveries",value=0),
                   numericInput("PC_whole","Proportion of Causal Variants (Optional)",value=NA),
                   numericInput("JJ_whole","Number of Variants (Optional)",value=NA)
                 ),
                 conditionalPanel(
                   condition = "input.SNPoption!='Whole Genome'",
                   numericInput("Alpha",HTML("&alpha;:"), value=0.0001),
                   selectInput("evoption", 
                               "Fixed single EV or Range of EV",
                               c("Single EV","Range of EV")),
                   conditionalPanel(
                     condition = "input.evoption == 'Single EV'",
                     numericInput("EV", "EV: Variance Explained(Percent)",
                                  value = 1)),
                   conditionalPanel(
                     condition = "input.evoption == 'Range of EV'",
                     sliderInput("EV_range","Variance Explained Range(Percent)",
                                 min=0.01,max=1,value = c(0.05,0.60),step=0.01)),
                   
                   conditionalPanel(
                     condition="input.SNPoption=='Single Gene'",
                     
                     selectInput(
                       "method", "Method",
                       c("SKAT","Calpha","Hotelling","Burden Test")
                     ),
                     
                     # Only show this panel if the plot type is a histogram
                     
                     conditionalPanel(
                       condition="input.method!='Burden Test'",
                       numericInput("PC","Proportion of Causal Variants (Optional)",value=NA)
                     ),   
                     conditionalPanel(
                       condition="input.method=='Burden Test'",
                       numericInput("PC_new","Proportion of Causal Variants",value=0.2),
                       numericInput("PRC","Proportion of Protective",value=0)
                     ),
                     
                     # textInput("nameEsseble","Specify Single Gene(Optional)"),
                     numericInput("JJ","Number of Variants (Optional)",value=NA)
                     
                   )
                 )
                 
                 
               ),
               conditionalPanel(
                 condition = "input.TypeofCalculation=='Sample Size Calculation'",
                 selectInput("SNPoption_s", 
                             "Type of power calculation",
                             c("Single Gene","Single SNP")),
                 
                 numericInput("PowerThreshold_s", "PowerThreshold", 
                              value=0.8),
    
                 conditionalPanel(
                   condition = "input.SNPoption_s!='Whole Genome'",
                   numericInput("Alpha_s",HTML("&alpha;:"), value=0.0001),
                   selectInput("evoption_s", 
                               "Fixed single EV or Range of EV",
                               c("Single EV","Range of EV")),
                   conditionalPanel(
                     condition = "input.evoption_s == 'Single EV'",
                     numericInput("EV_s", "EV: Variance Explained(Percent)",
                                  value = 1)),
                   conditionalPanel(
                     condition = "input.evoption_s == 'Range of EV'",
                     sliderInput("EV_range_s","Variance Explained Range(Percent)",
                                 min=0.01,max=1,value = c(0.05,0.60),step=0.01)),
                   
                   conditionalPanel(
                     condition="input.SNPoption=='Single Gene'",
                     
                     selectInput(
                       "method_s", "Method",
                       c("SKAT","Calpha","Hotelling","Burden Test")
                     ),
                     
                     # Only show this panel if the plot type is a histogram
                     
                     conditionalPanel(
                       condition="input.method_s!='Burden Test'",
                       numericInput("PC_s","Proportion of Causal Variants (Optional)",value=NA)
                     ),   
                     conditionalPanel(
                       condition="input.method_s=='Burden Test'",
                       numericInput("PC_new_s","Proportion of Causal Variants",value=0.2),
                       numericInput("PRC_s","Proportion of Protective",value=0)
                     ),
                     
                     # textInput("nameEsseble","Specify Single Gene(Optional)"),
                     numericInput("JJ_s","Number of Variants (Optional)",value=NA)
                     
                   )
                 )
                 
                 
                 
               ),
               
               
               actionButton("update","Update")
               
               
               
             ),
            
          
          
           
             # Show a table summarizing the values entered
             mainPanel(
               #tabsetPanel(
               #  tabPanel("Power Density",
               
               # p(HTML("Genetic Architecture I: &beta; independent of MAF"),style = "font-family: 'times'; font-si20pt"),
               #  p("Geneic Architecture II : EV independent of MAF",style = "font-family: 'times'; font-si20pt"), 
               # p(HTML("Genetic Architecture III: &beta;~log_10(MAF);"),style = "font-family: 'times'; font-si20pt"),
               conditionalPanel(
                 condition = "input.SNPoption=='Whole Genome'",
                 h5(" For whole genome power calculation, quick option runs genome-wide calculations  within 3 minutes and  provides rough estimates. Adequate option runs genome-wide calculations  within 6 minutes and provides more accurate estimates. Lastly, Complete option runs genome-wide calculations  within 15 minutes and  provides very accurate estimates.")),
               p("Senarario S1: MAF is independent of EV;",style = "font-family: 'times'; font-si20pt"),
               p("Senarario S2 : MAF is independent of genetic effects measured in the unit of per copy of an allele;",style = "font-family: 'times'; font-si20pt"), 
               p("Senarario S3: MAF is negatively correlated with genetic effect through the function ;",style = "font-family: 'times'; font-si20pt"),
               theme = "bootstrap1.css",
               busyIndicator(text = "Loading..", img = "busyIndicator/ajaxloaderq.gif", wait = 1),
              # DT::dataTableOutput("values"),
               plotOutput("plot")
               #   ),
               # tabPanel("Mean Power Distribution based on different EV",
               # h4("The result in this page will show up if you choose to specify a range of EV"),
               #    p(HTML("Genetic Architecture I: &beta; independent of MAF"),style = "font-family: 'times'; font-si20pt"),
               #   p("Geneic Architecture II : EV independent of MAF",style = "font-family: 'times'; font-si20pt"), 
               # p(HTML("Genetic Architecture III: &beta;~log_10(MAF);"),style = "font-family: 'times'; font-si20pt"),
               
               # dataTableOutput("values"),
               #dataTableOutput("values_r"),
               #dataTableOutput("values_r"),
               #plotOutput("plot_r")
               # )
               #)
             )
    ),
    tabPanel("Quantitative Trait",
             
             sidebarPanel( 
               selectInput("SNPoption2", 
                           "Type of power calculation",
                           c("Single Gene","Single SNP","Whole Genome")),
               
               numericInput("total2", "Total Sample Size:", 
                            value=10000),
               
               
               
               conditionalPanel(
                 condition = "input.SNPoption2=='Whole Genome'",
                 numericInput("Alpha_whole2",HTML("&alpha;:"), value=2.5e-06),
                 numericInput("GEV2", "GEV: Genome-wide Variance Explained(Percent)",
                              value = 20),
                 
                 selectInput(
                   "method_whole2", "Method",
                   c("SKAT","Calpha","Hotelling","Burden Test")
                 ),
                 selectInput("grid2",
                             "Computational Complexity(Optional)",
                             c("Quick","Adequate","Complete")),
                 numericInput("K2","Number of causal loci",value=1000),
                 numericInput("m2","Number of discoveries",value=0),
                 numericInput("PC_whole2","Proportion of Causal Variants (Optional)",value=NA),
                 numericInput("JJ_whole2","Number of Variants (Optional)",value=NA)
               ),
               conditionalPanel(
                 condition = "input.SNPoption2!='Whole Genome'",
                 numericInput("Alpha2",HTML("&alpha;:"), value=0.0001),
                 selectInput("evoption2", 
                             "Fixed single EV or Range of EV",
                             c("Single EV","Range of EV")),
                 conditionalPanel(
                   condition = "input.evoption2 == 'Single EV'",
                   numericInput("EV2", "EV: Variance Explained(Percent)",
                                value = 1)),
                 conditionalPanel(
                   condition = "input.evoption2 == 'Range of EV'",
                   sliderInput("EV_range2","Variance Explained Range(Percent)",
                               min=0.01,max=1,value = c(0.05,0.60),step=0.01)),
                 
                 conditionalPanel(
                   condition="input.SNPoption2=='Single Gene'",
                   
                   selectInput(
                     "method2", "Method",
                     c("SKAT","Calpha","Hotelling","Burden Test")
                   ),
                   
                   # Only show this panel if the plot type is a histogram
                   
                   conditionalPanel(
                     condition="input.method2!='Burden Test'",
                     numericInput("PC2","Proportion of Causal Variants (Optional)",value=NA)
                   ),   
                   conditionalPanel(
                     condition="input.method2=='Burden Test'",
                     numericInput("PC_new2","Proportion of Causal Variants",value=0.2),
                     numericInput("PRC2","Proportion of Protective",value=0)
                   ),
                   
                   # textInput("nameEsseble","Specify Single Gene(Optional)"),
                   numericInput("JJ2","Number of Variants (Optional)",value=NA)
                   
                 )
               ),
               
               
               actionButton("update2","Update")
             
             
             
          
             ),
             mainPanel(
             
               #tabsetPanel(
                # tabPanel("Power Density",
               conditionalPanel(
                 condition = "input.SNPoption2=='Whole Genome'",
                 h5(" For whole genome power calculation, quick option runs genome-wide calculations  within 3 minutes and  provides rough estimates. Adequate option runs genome-wide calculations  within 6 minutes and provides more accurate estimates. Lastly, Complete option runs genome-wide calculations  within 15 minutes and  provides very accurate estimates.")),
               p("Senarario S1: MAF is independent of EV;",style = "font-family: 'times'; font-si20pt"),
               p("Senarario S2 : MAF is independent of genetic effects measured in the unit of per copy of an allele;",style = "font-family: 'times'; font-si20pt"), 
               p("Senarario S3: MAF is negatively correlated with genetic effect through the function ;",style = "font-family: 'times'; font-si20pt"),
               theme = "bootstrap2.css",
               busyIndicator(text = "Loading..", img = "busyIndicator/ajaxloaderq.gif", wait = 1),
                          DT::dataTableOutput("values2"),
                          plotOutput("plot2")
               #  ),
                 #tabPanel("Mean Power Distribution based on different EV",
                          # h4("The result in this page will show up if you choose to specify a range of EV"),
                          # p(HTML("Genetic Architecture I: &beta; independent of MAF"),style = "font-family: 'times'; font-si20pt"),
                          # p("Geneic Architecture II : EV independent of MAF",style = "font-family: 'times'; font-si20pt"), 
                          # p(HTML("Genetic Architecture III: &beta;~log_10(MAF);"),style = "font-family: 'times'; font-si20pt"),
                          # dataTableOutput("values2_r"),
                          # plotOutput("plot2_r")
                 #)
              # )
             )
            
    ),

# tabPanel("Quantitative Trait Sample Size Calculation",
#          
#          sidebarPanel( 
#            selectInput("SNPoption2_s", 
#                        "Type of power calculation",
#                        c("Single Gene","Single SNP","Whole Genome")),
#            
#            numericInput("powerthreshold2_s", "Power Threshold:", 
#                         value=0.8),
#            
#            
#            
#            conditionalPanel(
#              condition = "input.SNPoption2_s=='Whole Genome'",
#              numericInput("Alpha_whole2_s",HTML("&alpha;:"), value=2.5e-06),
#              numericInput("GEV2_s", "GEV: Genome-wide Variance Explained(Percent)",
#                           value = 20),
#              
#              selectInput(
#                "method_whole2_s", "Method",
#                c("SKAT","Calpha","Hotelling","Burden Test")
#              ),
#              selectInput("grid2_s",
#                          "Computational Complexity(Optional)",
#                          c("Quick","Adequate","Complete")),
#              numericInput("K2_s","Number of causal loci",value=1000),
#              numericInput("m2_s","Number of discoveries",value=0),
#              numericInput("PC_whole2_s","Proportion of Causal Variants (Optional)",value=NA),
#              numericInput("JJ_whole2_s","Number of Variants (Optional)",value=NA)
#            ),
#            conditionalPanel(
#              condition = "input.SNPoption2!='Whole Genome'",
#              numericInput("Alpha2_s",HTML("&alpha;:"), value=0.0001),
#              selectInput("evoption2_s", 
#                          "Fixed single EV or Range of EV",
#                          c("Single EV","Range of EV")),
#              conditionalPanel(
#                condition = "input.evoption2_s == 'Single EV'",
#                numericInput("EV2_s", "EV: Variance Explained(Percent)",
#                             value = 1)),
#              conditionalPanel(
#                condition = "input.evoption2_s == 'Range of EV'",
#                sliderInput("EV_range2_s","Variance Explained Range(Percent)",
#                            min=0.01,max=1,value = c(0.05,0.60),step=0.01)),
#              
#              conditionalPanel(
#                condition="input.SNPoption2_s=='Single Gene'",
#                
#                selectInput(
#                  "method2_s", "Method",
#                  c("SKAT","Calpha","Hotelling","Burden Test")
#                ),
#                
#                # Only show this panel if the plot type is a histogram
#                
#                conditionalPanel(
#                  condition="input.method2_s!='Burden Test'",
#                  numericInput("PC2_s","Proportion of Causal Variants (Optional)",value=NA)
#                ),   
#                conditionalPanel(
#                  condition="input.method2_s=='Burden Test'",
#                  numericInput("PC_new2_s","Proportion of Causal Variants",value=0.2),
#                  numericInput("PRC2_s","Proportion of Protective",value=0)
#                ),
#                
#                # textInput("nameEsseble","Specify Single Gene(Optional)"),
#                numericInput("JJ2_s","Number of Variants (Optional)",value=NA)
#                
#              )
#            ),
#            
#            
#            actionButton("update2_s","Update")
#            
#            
#            
#            
#          ),
#          mainPanel(
#            
#            #tabsetPanel(
#            # tabPanel("Power Density",
#            conditionalPanel(
#              condition = "input.SNPoption2_s=='Whole Genome'",
#              h5(" For whole genome power calculation, quick option runs genome-wide calculations  within 3 minutes and  provides rough estimates. Adequate option runs genome-wide calculations  within 6 minutes and provides more accurate estimates. Lastly, Complete option runs genome-wide calculations  within 15 minutes and  provides very accurate estimates.")),
#            p("Senarario S1: MAF is independent of EV;",style = "font-family: 'times'; font-si20pt"),
#            p("Senarario S2 : MAF is independent of genetic effects measured in the unit of per copy of an allele;",style = "font-family: 'times'; font-si20pt"), 
#            p("Senarario S3: MAF is negatively correlated with genetic effect through the function ;",style = "font-family: 'times'; font-si20pt"),
#            theme = "bootstrap2.css",
#            busyIndicator(text = "Loading..", img = "busyIndicator/ajaxloaderq.gif", wait = 1),
#            DT::dataTableOutput("values2"),
#            plotOutput("plot2")
#            #  ),
#            #tabPanel("Mean Power Distribution based on different EV",
#            # h4("The result in this page will show up if you choose to specify a range of EV"),
#            # p(HTML("Genetic Architecture I: &beta; independent of MAF"),style = "font-family: 'times'; font-si20pt"),
#            # p("Geneic Architecture II : EV independent of MAF",style = "font-family: 'times'; font-si20pt"), 
#            # p(HTML("Genetic Architecture III: &beta;~log_10(MAF);"),style = "font-family: 'times'; font-si20pt"),
#            # dataTableOutput("values2_r"),
#            # plotOutput("plot2_r")
#            #)
#            # )
#          )
#          
# ),

    tabPanel("Citation",
             wellPanel(p("Derkach A., Zhang H. Chatterjee N. (2016)",a(href="http://biorxiv.org/content/early/2017/01/16/100891","Simplified Power Calculations for Aggregate-level Association Tests Provide Insights to Challenges for Rare Variant Association Studies.", "bioRxiv, 100891."))),            
    wellPanel(HTML('<a href="https://github.com/andrewhaoyu/ATTAP"><img src="img/githublogo.png" width=50 /> </a>')),
    h3("References"),
    p("Lek, M., Karczewski, K.J., Minikel, E.V., Samocha, K.E., Banks, E., Fennell, T., O'Donnell-Luria, A.H., Ware, J.S., Hill, A.J., Cummings, B.B., et al. (2016). Analysis of protein-coding genetic variation in 60,706 humans. Nature 536, 285-291.", style = "font-family: 'times'; font-si20pt"),
    p("Neale, B.M., Rivas, M.A., Voight, B.F., Altshuler, D., Devlin, B., Orho-Melander, M., Kathiresan, S., Purcell, S.M., Roeder, K., and Daly, M.J. (2011). Testing for an unusual distribution of rare variants. PLoS Genet 7, e1001322.
", style = "font-family: 'times'; font-si20pt"), 
    p("Derkach, A., Lawless, J.F., and Sun, L. (2013). Robust and powerful tests for rare variants using Fisher's method to combine evidence of association from two or more complementary tests. Genet Epidemiol 37, 110-121.", style = "font-family: 'times'; font-si20pt"), 
    p("Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M., and Lin, X. (2011). Rare-variant association testing for sequencing data with the sequence kernel association test. Am J Hum Genet 89, 82-93.", style = "font-family: 'times'; font-si20pt"), 
    p("Madsen, B.E., and Browning, S.R. (2009). A groupwise association test for rare mutations using a weighted sum statistic. PLoS Genet 5, e1000384.", style = "font-family: 'times'; font-si20pt"),
    p(" Li and Leal 2008, Am J Hum Genet. 2008 Sep;83(3):311-21. doi: 10.1016/j.ajhg.2008.06.024. Epub 2008 Aug 7.Methods for detecting associations with rare variants for common diseases: application to analysis of sequence data.
      ", style = "font-family: 'times'; font-si20pt"),
    p("  Neale et al. 2011, Derkach et al. 2013, Wu et al. 2011, Madsen and Browning 2009
      ", style = "font-family: 'times'; font-si20pt")
   
   
    
    )    
    
    
    
    
    
 
  )
)



