library(shiny)
library(shinythemes)
library(DT)


source('PowerCalc_Rare.r')
# Define UI for slider demo application
ui <- tagList(
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
               p("Power calculations are for \"most used\" gene based test called SKAT (http://www.hsph.harvard.edu/skat/).It provides average power. We average out all unknown parameters related to gene-phenotype architecture.Power calculations can be done also for linear + various weights, C-lapha, Hotelling, etc.", style = "font-family: 'times'; font-si20pt"),
               h3("Essential Parameters"),
               p("1.EV=% of variation in a trait explained by loci (gene),
                 2.Alpha = level of the test,
                 3.Number of cases,
                 4.Number of controls,
                 5.Total sample size
All you genes have ensembl names. So you can do power calculation for particular gene.Instead of name you can provide number of variants in a gene with MAF<0.01.", style = "font-family: 'times'; font-si20pt"),
               h3("Assumptions"),
               p("Right now, the power calculation is based on data from Exome Aggregation Consortium, so it is exome based calculations. It may underestimate number of variants in a gene (Based on write-up).
Power calculations can be done under two assumptions: 
                1) there is no relationship between MAF and % of variations explained by a variant 
                 2) there is no relationship between MAF and effect size (logOR). 
we have also developed method where one specify relationship between MAF and genetic effect beta", style = "font-family: 'times'; font-si20pt"),
               h3("Contributor"),
               p("AATAP is designed by Andriy Derkach, Haoyu Zhang and Nilanjan Chatterjee", style = "font-family: 'times'; font-si16pt")
             ,width = 15)
           )),
  tabPanel("Case Control",
  #  Application title
  #titlePanel("Aggregated Association test Approximate Power"),

    sidebarPanel(
      # Simple integer interval
      numericInput("ncases", "Number of Cases:", 
                value=5000),
      

      numericInput("ncont", "Number of Controls:", 
                  value = 5000),
      

      numericInput("EV", "EV: Variance Explained",
                  value = 0.01),
      

    numericInput("Alpha", "Alpha:", value=0.0019
                  ),
    textInput("nameEsseble","Specify One Gene Ensembl Name (optional)"),
    numericInput("JJ","Number of Variants (optional)",value=NULL),
    actionButton("update","Update")
    
      
    ),
    # Show a table summarizing the values entered
    mainPanel(
      verbatimTextOutput("explain"),
      dataTableOutput("values"),
      plotOutput("plot")
     )
),
  tabPanel("QT",
           sidebarPanel(
             # Simple integer interval
             numericInput("Total2", "Total Number of Subjects:", 
                          value=10000),
             
             
             numericInput("EV2", "EV: Variance Explained",
                          value = 0.01),
             
             
             numericInput("Alpha2", "Alpha:", value=0.0019
             ),
             numericInput("PC2","Percent of Causal (optional)",value=NULL),
             actionButton("update2","Update")
           ),
           mainPanel(
             verbatimTextOutput("explain2"),
             dataTableOutput("values2"),
             plotOutput("plot2")
           )
           
            
),
tabPanel("Citation")
)
)

server <- function(input, output) {
  
  # Reactive expression to compose a data frame containing all of
  # the values

  # Compose data frame


  #,input$nameEsseble,input$JJ



   data <- eventReactive(input$update,{
     if(input$nameEsseble=="" & is.na(input$JJ)){
       # round(get_Aprox(input$EV,
       #                input$Alpha,
       #                (input$ncases+input$ncont),
       #                input$ncases,input$ncont),
       #       3)
      get_Aprox(input$EV,
                input$Alpha,
                (input$ncases+input$ncont),
                input$ncases,input$ncont)
     }
     else if(input$nameEsseble=="" & is.na(input$JJ)!=1){
       # round(get_Aprox(input$EV,
       #                 input$Alpha,
       #                 (input$ncases+input$ncont),
       #                 input$ncases,
       #                 input$ncont,
       #                 JJ=input$JJ),
       #       3)
      get_Aprox(input$EV,
                       input$Alpha,
                       (input$ncases+input$ncont),
                       input$ncases,
                       input$ncont,
                       JJ=input$JJ)
         
     }
     else if(length(input$nameEsseble)!=0 & is.na(input$JJ) ){
       get_Aprox(input$EV,
                       input$Alpha,
                       (input$ncases+input$ncont),
                       input$ncases,
                       input$ncont,
                       nameEsseble=input$nameEsseble)

  
     }
     else
     {
       get_Aprox(input$EV,
                       input$Alpha,
                      (input$ncases+input$ncont),
                       input$ncases,input$ncont,
                       input$nameEsseble,
                       input$JJ)
     }
   })
   
   data2 <- eventReactive(input$update2,{
     if(is.na(input$JJ)){
       # round(get_Aprox(input$EV,
       #                input$Alpha,
       #                (input$ncases+input$ncont),
       #                input$ncases,input$ncont),
       #       3)
       get_Aprox(input$EV2,
                 input$Alpha2,
                 (input$Total2),
                 QT='QT')
     }
     else
     {
       get_Aprox(input$EV2,
                 input$Alpha2,
                 input$Total2,
                 QT='QT',
                 input$PC2
                 )
     }
   })

  # Show the values using an HTML table



  output$values <- DT::renderDataTable( 
    data.frame(
      Power_Distribution = c("Mean Power",
               "25% Quantile of Power",
               "Medium of Power",
               "75% Quantile of Power"
             ),
      Bin = c(
        data()[[1]][1],
        data()[[1]][2],
        data()[[1]][3],
        data()[[1]][4]
        ),
      Ein = c(
        data()[[1]][5],
        data()[[1]][6],
        data()[[1]][7],
        data()[[1]][8]
      ),
      Brel = c(
        data()[[1]][9],
        data()[[1]][10],
        data()[[1]][11],
        data()[[1]][12]
      )
    ),options=list(dom='t',autoWidth=T,columnDefs = list(list(width = '150px', targets = "_all"))),rownames=F
  )
  output$explain = renderText({("1.Bin: Beta independent of MAF
2.Ein : EV independent of MAF
3.Brel: Beta~log_10(MAF)"
                                )})
  output$plot = renderPlot(data()[[2]],height=300,width=600)
  
  output$explain2 = renderText({("1.Bin: Beta independent of MAF
2.Ein : EV independent of MAF
3.Brel: Beta~log_10(MAF)"
  )})
  
  output$values2 <- DT::renderDataTable( 
    data.frame(
      Power_Distribution = c("Mean Power",
                             "25% Quantile of Power",
                             "Medium of Power",
                             "75% Quantile of Power"
      ),
      Bin = c(
        data2()[[1]][1],
        data2()[[1]][2],
        data2()[[1]][3],
        data2()[[1]][4]
      ),
      Ein = c(
        data2()[[1]][5],
        data2()[[1]][6],
        data2()[[1]][7],
        data2()[[1]][8]
      ),
      Brel = c(
        data2()[[1]][9],
        data2()[[1]][10],
        data2()[[1]][11],
        data2()[[1]][12]
      )
    ),options=list(dom='t',autoWidth=T,columnDefs = list(list(width = '150px', targets = "_all"))),rownames=F
  )
  output$plot2 = renderPlot(data2()[[2]],height=300,width=600)
  
  
  
  
}

shinyApp(ui=ui,server=server)
