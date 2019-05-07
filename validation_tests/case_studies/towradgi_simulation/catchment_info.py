def create_catchment_list(simulation):

    from os.path import join
    
    channel_manning = simulation.args.channel_manning
    
    CatchmentList = [
        [join('Model', 'Bdy', 'Catchment.csv'), 100.0],
        [join('Model', 'Bdy', 'FineCatchment.csv'), 36.0],
        [join('Model', 'Bdy', 'CreekBanks.csv'), 8.0] 
        ]
    
    return CatchmentList


def create_manning_list(simulation):

    from os.path import join
    
    channel_manning = simulation.args.channel_manning
    
    ## IMPORTANT -- The ORDER in ManningList matters: When there is overlap,
    ##              priority regions at BOTTOM
    ## FIXME: This setup can be done with fewer lines of code!
    
    ManningList = [
       [ join('Model', 'Mannings', '1.csv'),0.04], #park
       [ join('Model', 'Mannings', '2.csv'),0.15],
       [ join('Model', 'Mannings', '3.csv'),0.15],
       [ join('Model', 'Mannings', '4.csv'),0.04],
       [ join('Model', 'Mannings', '5.csv'),0.15],
       [ join('Model', 'Mannings', '6.csv'),0.15],
       [ join('Model', 'Mannings', '7.csv'),0.15],
       [ join('Model', 'Mannings', '8.csv'),0.15],
       [ join('Model', 'Mannings', '9.csv'),0.04], #park
       [ join('Model', 'Mannings', '10.csv'), 0.15],
       [ join('Model', 'Mannings', '11.csv'), 0.15],
       [ join('Model', 'Mannings', '12.csv'), 0.15],
       [ join('Model', 'Mannings', '13.csv'), 0.04],
       [ join('Model', 'Mannings', '14.csv'), 0.15],
       [ join('Model', 'Mannings', '15.csv'), 0.15],
       [ join('Model', 'Mannings', '16.csv'), 0.15],
       [ join('Model', 'Mannings', '17.csv'), 0.15],
       [ join('Model', 'Mannings', '18.csv'), 0.045],
       [ join('Model', 'Mannings', '18a.csv'), 0.15],
       [ join('Model', 'Mannings', '18b.csv'), 0.15],
       [ join('Model', 'Mannings', '18c.csv'), 0.15],
       [ join('Model', 'Mannings', '18d.csv'), 0.15],
       [ join('Model', 'Mannings', '18e.csv'), 0.08], #cokeworks site
       [ join('Model', 'Mannings', '19.csv'), 0.15],
       [ join('Model', 'Mannings', '20.csv'), 0.15],
       [ join('Model', 'Mannings', '21.csv'), 0.15],
       [ join('Model', 'Mannings', '22.csv'), 0.15],
       [ join('Model', 'Mannings', '23.csv'), 0.15],
       [ join('Model', 'Mannings', '24.csv'), 0.05],
       [ join('Model', 'Mannings', '25.csv'), 0.15],
       [ join('Model', 'Mannings', '26.csv'), 0.15],
       [ join('Model', 'Mannings', '27.csv'), 0.15],
       [ join('Model', 'Mannings', '28.csv'), 0.15],
       [ join('Model', 'Mannings', '29.csv'), 0.15],
       [ join('Model', 'Mannings', '30.csv'), 0.15],
       [ join('Model', 'Mannings', '31.csv'), 0.15],
       [ join('Model', 'Mannings', '32.csv'), 0.15],
       [ join('Model', 'Mannings', '33.csv'), 0.15],
       [ join('Model', 'Mannings', '34.csv'), 0.15],
       [ join('Model', 'Mannings', '35.csv'), 0.15],
       [ join('Model', 'Mannings', '36.csv'), 0.05],
       [ join('Model', 'Mannings', '37.csv'), 0.15],
       [ join('Model', 'Mannings', '38.csv'), 0.15],
       [ join('Model', 'Mannings', '39.csv'), 0.15],
       [ join('Model', 'Mannings', '40.csv'), 0.15],
       [ join('Model', 'Mannings', '41.csv'), 0.15],
       [ join('Model', 'Mannings', '42.csv'), 0.15],
       [ join('Model', 'Mannings', '43.csv'), 0.15],
       [ join('Model', 'Mannings', '44.csv'), 0.15],
       [ join('Model', 'Mannings', '45.csv'), 0.15],
       [ join('Model', 'Mannings', '46.csv'), 0.15],
       [ join('Model', 'Mannings', '47.csv'), 0.15],
       [ join('Model', 'Mannings', '48.csv'), 0.15],
       [ join('Model', 'Mannings', '49.csv'), 0.15],
       [ join('Model', 'Mannings', '50.csv'), 0.15],
       [ join('Model', 'Mannings', '51.csv'), 0.15],
       [ join('Model', 'Mannings', '52.csv'), 0.15],
       [ join('Model', 'Mannings', '53.csv'), 0.15],
       [ join('Model', 'Mannings', '54.csv'), 0.15],
       [ join('Model', 'Mannings', '55.csv'), 0.15],
       [ join('Model', 'Mannings', '56.csv'), 0.15],
       [ join('Model', 'Mannings', '57.csv'), 0.15],
       [ join('Model', 'Mannings', '58.csv'), 0.15],
       [ join('Model', 'Mannings', '59.csv'), 0.08],
       [ join('Model', 'Mannings', '60.csv'), 0.15],
       [ join('Model', 'Mannings', '61.csv'), 0.08],
       [ join('Model', 'Mannings', '62.csv'), 0.15],
       [ join('Model', 'Mannings', '63.csv'), 0.08],
       [ join('Model', 'Mannings', '64.csv'), 0.15],
       [ join('Model', 'Mannings', '65.csv'), 0.15],
       [ join('Model', 'Mannings', '66.csv'), 0.15],
       [ join('Model', 'Mannings', '67.csv'), 0.15],
       [ join('Model', 'Mannings', '68.csv'), 0.15],
       [ join('Model', 'Mannings', '69.csv'), 0.15],
       [ join('Model', 'Mannings', '70.csv'), 0.15],
       [ join('Model', 'Mannings', '71.csv'), 0.05],
       [ join('Model', 'Mannings', '72.csv'), 0.15],
       [ join('Model', 'Mannings', '73.csv'), 0.15],
       [ join('Model', 'Mannings', '74.csv'), 0.15],
       [ join('Model', 'Mannings', '75.csv'), 0.15],
       [ join('Model', 'Mannings', '76.csv'), 0.15],
       [ join('Model', 'Mannings', '77.csv'), 0.07],
       [ join('Model', 'Mannings', '78.csv'), 0.15],
       [ join('Model', 'Mannings', '79.csv'), 0.15],
       [ join('Model', 'Mannings', '80.csv'), 0.15],
       [ join('Model', 'Mannings', '81.csv'), 0.15],
       [ join('Model', 'Mannings', '82.csv'), 0.15],
       [ join('Model', 'Mannings', '83.csv'), 0.15],
       [ join('Model', 'Mannings', '84.csv'), 0.15],
       [ join('Model', 'Mannings', '85.csv'), 0.15],
       [ join('Model', 'Mannings', '86.csv'), 0.15],
       [ join('Model', 'Mannings', 'Escarpement.csv'), 0.15],
       [ join('Model', 'Mannings', 'Railway.csv'), 0.04],
       [ join('Model', 'Creeks', 'creeks.csv'), channel_manning]
        ]

    return ManningList
