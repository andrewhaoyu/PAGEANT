library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)




source('PowerCalc_Rare.r')


# Define UI for slider demo application

function(input, output) {
  
  # Reactive expression to compose a data frame containing all of
  # the values
  
  # Compose data frame
  
  
  #,input$nameEsseble,input$JJ
  
  
  
  data <- eventReactive(input$update,{
    if(input$evoption=='No'){
    if((input$method)!='Burden Test'){
      if(is.na(input$JJ)&is.na(input$PC)){
        
        get_Aprox((input$EV)/100,
                  input$Alpha,
                  (input$ncases+input$ncont),
                  input$ncases,input$ncont,
                  TEST=input$method)
      }
      else if(is.na(input$JJ)&is.na(input$PC)!=1){
        
        get_Aprox((input$EV)/100,
                  input$Alpha,
                  (input$ncases+input$ncont),
                  input$ncases,
                  input$ncont,
                  PC=input$PC,
                  TEST=input$method)
        #return(data.frame(PC=input$PC,PRC=input$PRC))
        
      }
      else if(is.na(input$JJ)!=1&is.na(input$PC)){
        
        get_Aprox((input$EV)/100,
                  input$Alpha,
                  (input$ncases+input$ncont),
                  input$ncases,
                  input$ncont,
                  JJ=input$JJ,
                  TEST=input$method)
        
      }
      else{
        
        get_Aprox((input$EV)/100,
                  input$Alpha,
                  (input$ncases+input$ncont),
                  input$ncases,
                  input$ncont,
                  PC=input$PC,
                  TEST=input$method,
                  JJ=input$JJ)
        
      }
    
     
      
    }
    else if((input$method)=='Burden Test'){
      if(is.na(input$JJ)&is.na(input$PC_new)){
        
        get_Aprox((input$EV)/100,
                  input$Alpha,
                  (input$ncases+input$ncont),
                  input$ncases,input$ncont,
                  TEST=input$method,
                  PRC=input$PRC)
      }
      else if(is.na(input$JJ)&is.na(input$PC_new)!=1){
        
        get_Aprox((input$EV)/100,
                  input$Alpha,
                  (input$ncases+input$ncont),
                  input$ncases,
                  input$ncont,
                  TEST=input$method,
                  PC=input$PC_new,
                  PRC=input$PRC)
        
        #return(data.frame(PC=input$PC,PRC=input$PRC))
        
      }
      else if(is.na(input$JJ)!=1&is.na(input$PC_new)){
        
        get_Aprox((input$EV)/100,
                  input$Alpha,
                  (input$ncases+input$ncont),
                  input$ncases,
                  input$ncont,
                  JJ=input$JJ,
                  TEST=input$method,
                  PRC=input$PRC)
        
      }
      else{
        
        get_Aprox((input$EV)/100,
                  input$Alpha,
                  (input$ncases+input$ncont),
                  input$ncases,
                  input$ncont,
                  PC=input$PC_new,
                  TEST=input$method,
                  JJ=input$JJ,
                  PRC=input$PRC)
        
      }
     
    }
    
    }else{
    NULL
  }
    })

  
  data_r <- eventReactive(input$update,{
    if(input$evoption=='Yes'){
      if(input$method!='Burden Test'){
        
        if(is.na(input$JJ)&is.na(input$PC)){
          
          EV_new <- seq(input$EV_range[1],input$EV_range[2],(input$EV_range[2]-input$EV_range[1])/10)
          get_Aprox((EV_new)/100,
                    input$Alpha,
                    (input$ncases+input$ncont),
                    input$ncases,input$ncont,
                    TEST=input$method)
        }
        else if(is.na(input$JJ)&is.na(input$PC)!=1){
          
          EV_new <- seq(input$EV_range[1],input$EV_range[2],(input$EV_range[2]-input$EV_range[1])/10)
          get_Aprox(EV_new/100,
                    input$Alpha,
                    (input$ncases+input$ncont),
                    input$ncases,
                    input$ncont,
                    PC=input$PC,
                    TEST=input$method)
          
        }
        else if(is.na(input$JJ)!=1&is.na(input$PC)){
          
          EV_new <- seq(input$EV_range[1],input$EV_range[2],(input$EV_range[2]-input$EV_range[1])/10)
          get_Aprox(EV_new/100,
                    input$Alpha,
                    (input$ncases+input$ncont),
                    input$ncases,
                    input$ncont,
                    JJ=input$JJ,
                    TEST=input$method)
          
        }
        else{
          
          EV_new <- seq(input$EV_range[1],input$EV_range[2],(input$EV_range[2]-input$EV_range[1])/10)
          get_Aprox(EV_new/100,
                    input$Alpha,
                    (input$ncases+input$ncont),
                    input$ncases,
                    input$ncont,
                    PC=input$PC,
                    JJ=input$JJ,
                    TEST=input$method)
          
        }
       
        
      }
      else if(input$method=='Burden Test'){
        
        
        if(is.na(input$JJ)&is.na(input$PC_new)){
          
          EV_new <- seq(input$EV_range[1],input$EV_range[2],(input$EV_range[2]-input$EV_range[1])/10)
          get_Aprox((EV_new)/100,
                    input$Alpha,
                    (input$ncases+input$ncont),
                    input$ncases,input$ncont,
                    TEST=input$method,
                    PRC=input$PRC)
        }
        else if(is.na(input$JJ)&is.na(input$PC_new)!=1){
          
          EV_new <- seq(input$EV_range[1],input$EV_range[2],(input$EV_range[2]-input$EV_range[1])/10)
          get_Aprox(EV_new/100,
                    input$Alpha,
                    (input$ncases+input$ncont),
                    input$ncases,
                    input$ncont,
                    PC=input$PC_new,
                    TEST=input$method,
                    PRC=input$PRC)
          
        }
        else if(is.na(input$JJ)!=1&is.na(input$PC_new)!=1){
          
          EV_new <- seq(input$EV_range[1],input$EV_range[2],(input$EV_range[2]-input$EV_range[1])/10)
          get_Aprox(EV_new/100,
                    input$Alpha,
                    (input$ncases+input$ncont),
                    input$ncases,
                    input$ncont,
                    PC=input$PC_new,
                    JJ=input$JJ,
                    TEST=input$method,
                    PRC=input$PRC)
          
        }
       
     
        
      }
    }
    
  })
  
  data2 <- eventReactive(input$update2,{
    if(input$evoption2=='No'){
    if(input$method2!='Burden Test'){
      
      if(is.na(input$JJ2)&is.na(input$PC2)){
        
        get_Aprox((input$EV2)/100,
                  input$Alpha2,
                  (input$total2),
                  QT='QT',
                  TEST=input$method2)
      }
      else if(is.na(input$JJ2)&is.na(input$PC2)!=1){
        
        get_Aprox(input$EV2/100,
                  input$Alpha2,
                  (input$total2),
                  PC=input$PC2,
                  QT='QT',
                  TEST=input$method2)
        
      }
      else if(is.na(input$JJ2)!=1&is.na(input$PC2)){
        
        get_Aprox(input$EV2/100,
                  input$Alpha2,
                  (input$total2),
                  JJ=input$JJ2,
                  QT='QT',
                  TEST=input$method2)
        
      }
      else if(is.na(input$JJ2)!=1&is.na(input$PC2)!=1){
        
        get_Aprox(input$EV2/100,
                  input$Alpha2,
                  (input$total),
                  PC=input$PC2,
                  JJ=input$JJ2,
                  QT='QT',
                  TEST=input$method2)
        
      }
     
      
    }
    else if(input$method2=='Burden Test'){
      
      
      if(is.na(input$JJ2)&is.na(input$PC_new2)){
        
        get_Aprox((input$EV2)/100,
                  input$Alpha2,
                  (input$total2),
                  QT='QT',
                  TEST=input$method2,
                  PRC=input$PRC2)
      }
      else if( is.na(input$JJ2)&is.na(input$PC_new2)!=1){
        
        get_Aprox(input$EV2/100,
                  input$Alpha2,
                  (input$total2),
                  PC=input$PC_new2,
                  QT='QT',
                  TEST=input$method2,
                  PRC=input$PRC2)
        
      }
      else if(is.na(input$JJ2)!=1&is.na(input$PC_new2)){
        
        get_Aprox(input$EV2/100,
                  input$Alpha2,
                  (input$total2),
                  JJ=input$JJ2,
                  QT='QT',
                  TEST=input$method2,
                  PRC=input$PRC2)
        
      }
      else if(is.na(input$JJ2)!=1&is.na(input$PC_new2)!=1){
        
        get_Aprox(input$EV2/100,
                  input$Alpha2,
                  (input$total),
                  PC=input$PC2_new2,
                  JJ=input$JJ2,
                  QT='QT',
                  TEST=input$method2,
                  PRC=input$PRC2)
        
      }
      else if( is.na(input$JJ2)&is.na(input$PC_new2)){
        
        get_Aprox(input$EV2/100,
                  input$Alpha2,
                  (input$total2),
                  nameEsseble = input$nameEsseble2,
                  QT='QT',
                  TEST=input$method2,
                  PRC=input$PRC2)
        
      }
      
      
    }
    }else{
    
  }
  
})
  
  
  data2_r <- eventReactive(input$update2,{
    if(input$evoption2=='Yes'){
      if(input$method2!='Burden Test'){
        
        if(is.na(input$JJ2)&is.na(input$PC2)){
          
          EV_new <- seq(input$EV_range2[1],input$EV_range2[2],(input$EV_range2[2]-input$EV_range2[1])/10)
          get_Aprox(EV_new/100,
                    input$Alpha2,
                    (input$total2),
                    QT='QT',
                    TEST=input$method2)
        }
        else if(is.na(input$JJ2)&is.na(input$PC2)!=1){
          EV_new <- seq(input$EV_range2[1],input$EV_range2[2],(input$EV_range2[2]-input$EV_range2[1])/10)
          get_Aprox(EV_new/100,
                    input$Alpha2,
                    (input$total2),
                    PC=input$PC2,
                    QT='QT',
                    TEST=input$method2)
          
        }
        else if(is.na(input$JJ2)!=1&is.na(input$PC2)){
          EV_new <- seq(input$EV_range2[1],input$EV_range2[2],(input$EV_range2[2]-input$EV_range2[1])/10)
          get_Aprox(EV_new/100,
                    input$Alpha2,
                    (input$total2),
                    JJ=input$JJ2,
                    QT='QT',
                    TEST=input$method2)
          
        }
        else if(is.na(input$JJ2)!=1&is.na(input$PC2)!=1){
          EV_new <- seq(input$EV_range2[1],input$EV_range2[2],(input$EV_range2[2]-input$EV_range2[1])/10)
          get_Aprox(EV_new/100,
                    input$Alpha2,
                    (input$total),
                    PC=input$PC2,
                    JJ=input$JJ2,
                    QT='QT',
                    TEST=input$method2)
          
        }
       
        
      }
      else if(input$method2=='Burden Test'){
        
        
        if(is.na(input$JJ2)&is.na(input$PC_new2)){
          
          EV_new <- seq(input$EV_range2[1],input$EV_range2[2],(input$EV_range2[2]-input$EV_range2[1])/10)
          get_Aprox(EV_new/100,
                    input$Alpha2,
                    (input$total2),
                    QT='QT',
                    TEST=input$method2,
                    PRC=input$PRC2)
        }
        else if(is.na(input$JJ2)&is.na(input$PC_new2)!=1){
          EV_new <- seq(input$EV_range2[1],input$EV_range2[2],(input$EV_range2[2]-input$EV_range2[1])/10)
          get_Aprox(EV_new/100,
                    input$Alpha2,
                    (input$total2),
                    PC=input$PC_new2,
                    QT='QT',
                    TEST=input$method2,
                    PRC=input$PRC2)
          
        }
        else if(is.na(input$JJ2)!=1&is.na(input$PC_new2)){
          EV_new <- seq(input$EV_range2[1],input$EV_range2[2],(input$EV_range2[2]-input$EV_range2[1])/10)
          get_Aprox(EV_new/100,
                    input$Alpha2,
                    (input$total2),
                    JJ=input$JJ2,
                    QT='QT',
                    TEST=input$method2,
                    PRC=input$PRC2)
          
        }
        else if(is.na(input$JJ2)!=1&is.na(input$PC_new2)!=1){
          EV_new <- seq(input$EV_range2[1],input$EV_range2[2],(input$EV_range2[2]-input$EV_range2[1])/10)
          get_Aprox(EV_new/100,
                    input$Alpha2,
                    (input$total),
                    PC=input$PC_new2,
                    JJ=input$JJ2,
                    QT='QT',
                    TEST=input$method2,
                    PRC=input$PRC2)
          
        }
       
        
        
      }
    }
  })
  
  
  # Show the values using an HTML table
  
  
    output$values <- DT::renderDataTable( 
      withProgress(message="We are computing",{ 
      data()[[1]]
         })  
      
      ,options=list(dom='t',autoWidth=T,columnDefs = list(list(width = '150px', targets = "_all"))),rownames=F
    )

    output$plot = renderPlot(data()[[2]],height=300,width=600)
    
    
    output$values_r <- DT::renderDataTable( 
      withProgress(message="We are based on different EV range",{ 
        data_r()[[1]]
         })  
      
      ,options=list(dom='t',autoWidth=T,columnDefs = list(list(width = '150px', targets = "_all"))),rownames=F
    )
    
    
    #output$explain2 = renderText({("Genetic Architecture I: &beta; independent of MAF
#Geneic Architecture II : EV independent of MAF
#Genetic Architecture III: &beta;~log_10(MAF)")
    #})
   
   output$plot_r <- renderPlot(data_r()[[2]],height = 300,width = 600)
   
   
     
   output$values2 <- DT::renderDataTable( 
     withProgress(message="We are computing",{ 
   data2()[[1]]
        })  
     
     ,options=list(dom='t',autoWidth=T,columnDefs = list(list(width = '150px', targets = "_all"))),rownames=F
   )
   output$plot2 = renderPlot(data2()[[2]],height=300,width=600)
    
   output$values2_r <- DT::renderDataTable( 
     withProgress(message="We are computing based on different EV range",{ 
       data2_r()[[1]]
     })  
     ,options=list(dom='t',autoWidth=T,columnDefs = list(list(width = '150px', targets = "_all"))),rownames=F
   )
   output$plot2_r <- renderPlot(data2_r()[[2]],height = 300,width = 600)
    
 
  
  
  
  }