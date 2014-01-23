
shinyServer(function(input, output) {
  
  ## Main reactive function
  mod <- reactive({
    epiDCM(
        type = input$modtype, 
        groups = 1, 
        s.num = input$s.num, 
        i.num = input$i.num, 
        r.num = input$r.num, 
        trans.rate = input$trans.rate, 
        act.rate = input$act.rate, 
        rec.rate = input$rec.rate, 
        b.rate = input$b.rate, 
        ds.rate = input$ds.rate, 
        di.rate = input$di.rate, 
        dr.rate = input$dr.rate, 
        nsteps = input$nsteps,
        dt = input$dt
      )
  })
  
  ## Main Plot tab
  output$MainPlot <- renderPlot({
    par(mar=c(3.5,3.5,1.2,1), mgp=c(2.1,1,0))
    if (input$compsel == 'State Prevalence')
      plot(mod(), leg.cex=1.1, alpha=input$alpha, lwd=3.5)
    if (input$compsel == 'State Size')
      plot(mod(), popfrac=F, leg.cex=1.1, alpha=input$alpha, lwd=3.5)
    if (input$compsel == 'Incidence')
      plot(mod(), y='si.flow', popfrac=F, leg.cex=1.1, alpha=input$alpha, lwd=3.5)
  })
  output$dlMainPlot <- downloadHandler(
    filename = 'MainPlot.pdf',
    content = function(file) {
      pdf(file=file, h=6, w=10)
      par(mar=c(3.5,3.5,1.2,1), mgp=c(2.1,1,0))
      if (input$compsel == 'State Prevalence')
        plot(mod(), leg.cex=1.1, alpha=input$alpha, lwd=3.5)
      if (input$compsel == 'State Size')
        plot(mod(), popfrac=F, leg.cex=1.1, alpha=input$alpha, lwd=3.5)
      if (input$compsel == 'Incidence')
        plot(mod(), y='si.flow', popfrac=F, leg.cex=1.1, alpha=input$alpha, lwd=3.5)
      dev.off()
    }
  )

  ## Summary and Compartment plot tab
  # Outfrom from summary
  output$outSummary <- renderPrint({
    summary(mod(), time=input$summTs, digits=input$summDig)
  })

  # Comp.plot
  output$CompPlot <- renderPlot({
    comp.plot(mod(), time=input$summTs, digits=input$summDig)
  })
  
  # Download for comp.plot
  output$dlCompPlot <- downloadHandler(
    filename = 'CompPlot.pdf',
    content = function(file) {
      pdf(file=file, h=6, w=10)
      comp.plot(mod(), time=input$summTs, digits=input$summDig)
      dev.off()
    }
  )
  
  
  ## Data tab
  # Stats tab
  output$dlData <- downloadHandler(
    filename = 'ModelData.csv',
    content = function(file) {
      write.csv(as.data.frame(mod()), file)
    }
  )
  output$outData <- renderTable({
    as.data.frame(mod())
  })
  

})