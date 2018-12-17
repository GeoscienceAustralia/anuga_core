
import numpy

from anuga import Boyd_box_operator
from os.path import join



def setup_structures(simulation):
    """
    Setup Culvert Data
    """
    
    domain  = simulation.domain
    args = simulation.args
    channel_manning = args.channel_manning
    
    """
    #FIXME: Include these again
    
    #------------------------------------------------------------------------------
    # ENTER CULVERT DATA
    #------------------------------------------------------------------------------
    # BrookerStCulvert branch 2
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[305772.982,6193988.557] , [305772.378,6193987.823]])
    el1 = numpy.array([[305794.592,6193983.907] , [305793.988,6193983.173]])
    
    ## Adjust el0, el1
    #elOffset=0.
    #el0M=0.5*(el0[0,:]+el0[1,:]) ; el1M=0.5*(el1[0,:]+el1[1,:]); n0=el0M-el1M; n0=n0/((n0*n0).sum())**0.5;
    #el0 = el0    
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=0.9,
                                exchange_lines=[el0, el1],
                                height=0.9,
                                apron=3.0,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '1H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)    
    
    # MeadowStCulvert branch 2
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[305886.333,6193929.052] , [305883.172,6193922.986]])
    el1 = numpy.array([[305906.553,6193910.461] , [305903.393,6193904.395]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=5.4,
                                exchange_lines=[el0, el1],
                                height=0.6,
                                apron=3.0,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '2H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)    
    
    # WilliamsStCulvert branch 2
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[305945.955,6193836.293] , [305945.125,6193835.387]])
    el1 = numpy.array([[306040.565,6193827.573] , [306039.735,6193826.667]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=1.2,
                                exchange_lines=[el0, el1],
                                height=1.2,
                                apron=3.0,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '3H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)     
    
    # MeadowStCulvert Towradgi branch
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[305812.113,6193591.972] , [305809.390,6193588.820]])
    el1 = numpy.array([[305834.913,6193588.382] , [305832.190,6193585.230]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=4.0,
                                exchange_lines=[el0, el1],
                                height=2.2,
                                apron=3.0,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '4H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)  
    
    # CollinsStCulverts tarra 5 branch
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[306330.608,6194817.116] , [306320.768,6194805.884]])
    el1 = numpy.array([[306369.483,6194811.616] , [306359.643,6194800.384]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=3.6,
                                exchange_lines=[el0, el1],
                                height=0.93,
                                apron=3.0,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '9H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)                                     
                                        
    # Norther Distributor Culverts branch 5
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[306956.242,6194465.589] , [306950.446,6194457.411]])
    el1 = numpy.array([[307003.711,6194446.089] , [306997.916,6194437.911]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=9.0,
                                exchange_lines=[el0, el1],
                                height=0.85,
                                apron=3.0,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '10H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)                                      
                                        
    #CokeWorksCulverts branch 5
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[307142.161,6194181.3065] , [307138.519,6194174.394]])
    el1 = numpy.array([[307160.521,6194164.8165] , [307156.879,6194157.904]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=6.2,
                                exchange_lines=[el0, el1],
                                height=3.0,
                                apron=3.1,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '11H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)                                      
                                        
    #Northern Distributor Branch 6 Culverts
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[306950.758,6193454.717] , [306947.804,6193453.283]])
    el1 = numpy.array([[306988.633,6193474.217] , [306985.679,6193472.783]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=3.6,
                                exchange_lines=[el0, el1],
                                height=1.2,
                                apron=3.1,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '12H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)                                      
                                        
    #Railway Branch 6 Culverts
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[307139.134,6193474.458] , [307138.492,6193473.542]])
    el1 = numpy.array([[307150.884,6193469.458] , [307150.242,6193468.542]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=3.6,
                                exchange_lines=[el0, el1],
                                height=1.2,
                                apron=3.1,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '13H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False) 
    
    #Colgong Branch6 Culverts
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[307200.610,6193476.765] , [307199.140,6193475.235]])
    el1 = numpy.array([[307224.610,6193475.765] , [307223.140,6193474.235]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=2.1,
                                exchange_lines=[el0, el1],
                                height=1.0,
                                apron=3.1,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '14H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)   
                                        
    #Basin Outlet Branch 3 Culverts
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[305629.639,6194408.883] , [305626.521,6194400.457]])
    el1 = numpy.array([[305665.889,6194347.183] , [305662.771,6194338.757]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=6.0,
                                exchange_lines=[el0, el1],
                                height=0.86,
                                apron=3.1,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '15H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)                                      
                                        
    #BellambiRd Branch 3 Culverts
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[305777.182,6194305.377] , [305776.444,6194304.623]])
    el1 = numpy.array([[305873.807,6194303.377] , [305873.069,6194302.623]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=1.05,
                                exchange_lines=[el0, el1],
                                height=1.0,
                                apron=3.1,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '16H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)    
    
    #MeadowSt Branch 3 Culverts
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[305914.649,6194322.375] , [305913.477,6194321.625]])
    el1 = numpy.array([[305950.711,6194335.375] , [305949.539,6194334.625]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=1.5,
                                exchange_lines=[el0, el1],
                                height=1.0,
                                apron=3.1,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '17H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)     
    
    #13MeadowSt Branch 3 Culverts
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[305911.280,6194359.203] , [305910.260,6194358.017]])
    el1 = numpy.array([[305946.090,6194353.573] , [305945.070,6194352.387]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=1.5,
                                exchange_lines=[el0, el1],
                                height=1.0,
                                apron=3.1,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '18H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)     
    
    #41AngelStBranch3Culverts
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[306196.779,6194028.193] , [306192.221,6194010.807]])
    el1 = numpy.array([[306200.154,6194018.693] , [306195.596,6194001.307]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=10.0,
                                exchange_lines=[el0, el1],
                                height=0.35,
                                apron=3.1,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '19H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)        
    
    #CarrollSt Branch 7 Culverts
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[308002.045,6193820.163] , [308001.215,6193819.197]])
    el1 = numpy.array([[308021.965,6193816.883] , [308021.135,6193815.917]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=1.22,
                                exchange_lines=[el0, el1],
                                height=0.3,
                                apron=3.1,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '20H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)           
    
    #ParkerRd Branch 7 Culverts
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[308105.832,6193803.622] , [308103.648,6193801.118]])
    el1 = numpy.array([[308126.782,6193800.552] , [308124.598,6193798.048]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=3.2,
                                exchange_lines=[el0, el1],
                                height=0.3,
                                apron=3.1,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '21H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)     
    
    #LakePde Branch 7 Culverts
    losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
    el0 = numpy.array([[308251.257,6193614.658] , [308248.343,6193618.]])
    el1 = numpy.array([[308232.,6193593.] , [308225.,6193596.]])
        
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=3.0,
                                exchange_lines=[el0, el1],
                                height=0.75,
                                apron=3.1,
                                enquiry_gap=10.0,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=0.013,
                                logging=False,
                                label=join('Runs', '22H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)                                      
    """                                   

    
    
    smoothTS=30. # Smoothing timescale for bridges
    # PrincesHwyBridge Towradgi branch
    losses = {'inlet':0.0, 'outlet':0.0, 'bend':0.0, 'grate':0.0, 'pier': 1.0, 'other': 0.0}
    el0 = numpy.array([[306608.,6193703.] , [306607.,6193700.0]]) 
    el1 = numpy.array([[306624.,6193692.7] , [306623.,6193688.]])
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=12.0,
                                exchange_lines=[el0, el1],
                                height=3.0,
                                apron=0.0,
                                enquiry_gap=10.0,
                                smoothing_timescale=smoothTS,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=channel_manning,
                                logging=False,
                                label=join('Runs', '5H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)  
    
    # PioneerRdBridge Towradgi branch
    losses = {'inlet':0.0, 'outlet':0.0, 'bend':0.0, 'grate':0.0, 'pier': 1.0, 'other': 0.0}
    el0 = numpy.array([[307623.,6193610.] , [307622.,6193607.]])
    el1 = numpy.array([[307610.,6193619.] , [307609., 6193616.]])  
    enq0 = numpy.array([[ 307637., 6193588. ]])  
    enq1 = numpy.array([[ 307596., 6193623. ]])  
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=20.0,
                                exchange_lines=[el0, el1],
                                #enquiry_points=[enq0, enq1], # This seemed to make stability worse
                                height=3.5,
                                apron=0.0,
                                enquiry_gap=10.0,
                                smoothing_timescale=smoothTS,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=channel_manning,
                                logging=False,
                                label=join('Runs', '8H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)                           
    
    # NorthernDistributorBridge Towradgi branch
    losses = {'inlet':0.0, 'outlet':0.0, 'bend':0.0, 'grate':0.0, 'pier': 1.0, 'other': 0.0}
    el0 = numpy.array([[306985.,6193749.] , [306985.,6193736.]])
    el1 = numpy.array([[306950.,6193745.] , [306950.,6193732.]])
    
    #enq0 = numpy.array([[ 306996., 6193750.]])
    #enq1 = numpy.array([[ 306931., 6193738.]])
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=45.0,
                                exchange_lines=[el0, el1],
                                #enquiry_points=[enq0, enq1],
                                height=6.0,
                                apron=0.0,
                                enquiry_gap=10., #None, #10.0,
                                smoothing_timescale=smoothTS,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=channel_manning,
                                logging=False,
                                label=join('Runs', '6H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False)    
    
    # RailwayBridge Towradgi branch ? (THIS IS THE BRIDGE THAT BLOCKEDAND DAM BREAK OCCURRED) ?
    losses = {'inlet':0.0, 'outlet':0.0, 'bend':0.0, 'grate':0.0, 'pier': 1.0, 'other': 0.0}
    el0 = numpy.array([[307236.,6193737.] , [307235.,6193733.]])
    el1 = numpy.array([[307223.,6193738.] , [307222.,6193734.]]) 
    culvert = Boyd_box_operator(domain,
                                losses=losses,
                                width=20.0,
                                exchange_lines=[el0, el1],
                                height=8.0,
                                apron=0.0,
                                enquiry_gap=20.0,
                                smoothing_timescale=smoothTS,
                                use_momentum_jet=True,
                                use_velocity_head=True,
                                manning=channel_manning,
                                logging=False,
                                label=join('Runs', '7H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                verbose=False) 
