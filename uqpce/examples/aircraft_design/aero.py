import numpy as np
import openmdao.api as om
#hi
#dtermined inputs
#planform area (S)
# Aspect Ratio (AR)
# Crusie Speed (V_cruise)

#coupled inputs
# m_total

#uncertain inputs (aleatory, I think)
# delta_CD0, delta_ks, delta_e

#outputs 
# LD, CL, CD

from fixed import parameters


from scipy.special import erfinv, erf
import matplotlib.pyplot as plt

class AeroDiscipline(om.ExplicitComponent):

    def initialize(self):
        self.options.declare('vec_size', types=int)

    def setup(self):
        n = self.options['vec_size']
       
        self.add_input('g', val=parameters['g'], units="m/s**2" )
        self.add_input('rho', val=parameters['rho'], units="kg/m**3")
        self.add_input('C_D0_base', val=parameters['CD0_base'], units=None)
        self.add_input('S_0', val=parameters['S_naught'], units="m**2" )
        self.add_input('ks_base',val=parameters['ks_base'], units="1/m**2")
        self.add_input('e_base', val=parameters['e_oswald_base'], units=None)
        
        self.add_input('S', val=parameters['S'], units="m**2")
        self.add_input('V', val=parameters['V'], units="m/s")
        self.add_input('AR', val=parameters['AR'],units=None)
    
        self.add_input('m_total',val=parameters['m_total'],units="kg",
                       shape=(n,))
        self.add_input('delta_CD0',val=1.0,units=None,
                       shape=(n,))
        self.add_input('delta_ks',val=1.0,units=None,
                       shape=(n,))
        self.add_input('delta_e',val=1.0,units=None,
                       shape=(n,))
    
        self.add_output('CL',0.0,shape=(n,),units=None)
        self.add_output('CD',0.0,shape=(n,),units=None)
        self.add_output('LD',0.0,shape=(n,),units=None)
        self.add_output('WL',0.0,shape=(n,),units=None)

    
    def setup_partials(self):
        n = self.options['vec_size']
        arange = np.arange(n)
    
    #Sensitivities-start~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        self.declare_partials(of="CL",wrt="V",method="exact")
        self.declare_partials(of="CL",wrt="S",method="exact")
        #self.declare_partials(of="CL",wrt="AR",method="exact")
        #des variables^
        self.declare_partials(of="CL",wrt="m_total",method="exact", rows=arange, cols=arange)
        self.declare_partials(of="CL",wrt="rho",method="exact")
        self.declare_partials(of="CL",wrt="g",method="exact")
        #all the rest of CL wrt other inputs are zero by default, so not needed
        #other partials just in case^

        self.declare_partials(of="CD",wrt="V",method="exact")
        self.declare_partials(of="CD",wrt="S",method="exact")
        self.declare_partials(of="CD",wrt="AR",method="exact")
        self.declare_partials(of="CD",wrt="m_total",method="exact", rows=arange, cols=arange)
        self.declare_partials(of="CD",wrt="rho",method="exact")
        self.declare_partials(of="CD",wrt="g",method="exact")
        self.declare_partials(of="CD",wrt="C_D0_base",method="exact")
        self.declare_partials(of="CD",wrt="S_0",method="exact")
        self.declare_partials(of="CD",wrt="e_base",method="exact")
        self.declare_partials(of="CD",wrt="delta_CD0",method="exact", rows=arange, cols=arange)
        self.declare_partials(of="CD",wrt="delta_e",method="exact", rows=arange, cols=arange)
        self.declare_partials(of="CD",wrt="ks_base",method="exact")
        self.declare_partials(of="CD",wrt="delta_ks",method="exact", rows=arange, cols=arange)


        self.declare_partials(of="LD",wrt="V",method="exact")
        self.declare_partials(of="LD",wrt="S",method="exact")
        self.declare_partials(of="LD",wrt="AR",method="exact")
        self.declare_partials(of="LD",wrt="m_total",method="exact", rows=arange, cols=arange)
        self.declare_partials(of="LD",wrt="rho",method="exact")
        self.declare_partials(of="LD",wrt="g",method="exact")
        self.declare_partials(of="LD",wrt="C_D0_base",method="exact")
        self.declare_partials(of="LD",wrt="S_0",method="exact")
        self.declare_partials(of="LD",wrt="e_base",method="exact")
        self.declare_partials(of="LD",wrt="delta_CD0",method="exact", rows=arange, cols=arange)
        self.declare_partials(of="LD",wrt="delta_e",method="exact", rows=arange, cols=arange)
        self.declare_partials(of="LD",wrt="ks_base",method="exact")
        self.declare_partials(of="LD",wrt="delta_ks",method="exact", rows=arange, cols=arange)

        self.declare_partials(of="WL",wrt="S",method="exact", rows=arange, cols=arange)
        self.declare_partials(of="WL",wrt="m_total",method="exact", rows=arange, cols=arange)
        self.declare_partials(of="WL",wrt="g",method="exact", rows=arange, cols=arange)


    #Sensitivities-end~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #passes input member inherited from om.Exp for reading and
    #outputs memeber struct/map thing whatever python calls it for writing
    def compute(self,inputs,outputs):
        g = inputs['g']
        rho = inputs['rho']
        C_D0_base = inputs['C_D0_base']
        S_0 = inputs['S_0']
        ks_base = inputs['ks_base']
        e_base = inputs['e_base']
        delta_CD0 = inputs['delta_CD0']
        delta_ks = inputs['delta_ks']
        delta_e = inputs['delta_e']


        #do this \/ double equal thingy to reuse output when needed, this synatx pattern might be useful
        #in compute partials function for chain rule stuff
        outputs['CL'] = CL = (inputs['m_total']*g) / ((1.0/2.0)*rho*(inputs['V']**2)*inputs['S'])
        C_D0 = C_D0_base*delta_CD0 + ks_base*delta_ks*(inputs['S']-S_0)     
        outputs['CD'] = CD = C_D0 + (CL**2) / (np.pi*inputs['AR']*e_base*delta_e)
        outputs['LD'] = CL/CD
        outputs['WL'] = (inputs['m_total']*g) / inputs['S']

    def compute_partials(self, inputs, partials): #I presume inputs and partials are inherited memebers of
        g = inputs['g']
        rho = inputs['rho']
        C_D0_base = inputs['C_D0_base']
        S_0 = inputs['S_0']
        ks_base = inputs['ks_base']
        e_base = inputs['e_base']
        delta_CD0 = inputs['delta_CD0']
        delta_ks = inputs['delta_ks']
        delta_e = inputs['delta_e']

        CL = (inputs['m_total']*g) / ((1.0/2.0)*rho*(inputs['V']**2)*inputs['S'])         
        C_D0 = C_D0_base*delta_CD0 + ks_base*delta_ks*(inputs['S']-S_0) 
        CD = C_D0 + (CL**2) / (np.pi*inputs['AR']*e_base*delta_e)
                                                  
        partials['CL','V'] = dCLdV = -2*CL*(1.0/inputs['V'])
        partials['CL','S'] = dCLdS = -1*CL*(1.0/inputs['S'])
        #partials['CL','AR'] = dCLdAR = 0 #fixed to assume S and AR as independent. span is always 
        #computed from these inputs
        dCLdAR = 0
        partials['CL','m_total'] = dCLdm = CL/inputs['m_total']
        partials['CL','rho'] = dCLdrho = -CL/rho
        partials['CL','g'] = dCLdg = CL/g
    
        #ugliness helpers
        dCD_0dV = 0.0
        dCD_0dS = ks_base*delta_ks
        #b_squared = inputs['AR']*inputs['S']
        dSdAR = 0.0 #fixed to assume S and AR as independent. span is always 
        #computed from these inputs
        dCD_0dAR = dCD_0dS*dSdAR
        dARdS = 0.0 #fixed to assume S and AR as independent. span is always 
        #computed from these inputs
        dCD_0dm = 0.0
        dCD_0drho = 0.0 
        dCD_0dg = 0.0
        #dARdg = 0.0
        dCD_0dCDbase = delta_CD0
        dCD_0dS0 = -ks_base*delta_ks
        dCD_0debase = 0.0
        dCD_0ddeltaCD0 = C_D0_base
        dCD_0ddeltae = 0


        #product rule/quotient rule or whatever u wanna call it helpers
        product_rule_V = 2*CL*dCLdV*(1/inputs['AR']) #+ (CL**2)*(0)
        product_rule_S = 2*CL*dCLdS*(1/inputs['AR']) - (CL**2)*(1.0/(inputs['AR']**2))*dARdS
        product_rule_AR = 2*CL*dCLdAR*(1/inputs['AR']) - (CL**2)*(1.0/(inputs['AR']**2))*(1.0)
        product_rule_m = 2*CL*dCLdm*(1/inputs['AR'])
        product_rule_rho = 2*CL*dCLdrho*(1/inputs['AR'])
        product_rule_g = 2*CL*dCLdg*(1/inputs['AR'])
        
        partials['CD','V'] = dCDdV = dCD_0dV + (1/(np.pi*e_base*delta_e))*(product_rule_V)
        partials['CD','S'] = dCDdS =  dCD_0dS + (1/(np.pi*e_base*delta_e))*(product_rule_S)
        partials['CD','AR'] = dCDdAR = dCD_0dAR +  (1/(np.pi*e_base*delta_e))*(product_rule_AR)
        partials['CD','m_total'] = dCDdm =  dCD_0dm + (1/(np.pi*e_base*delta_e))*(product_rule_m)
        partials['CD','rho'] = dCDdrho =  dCD_0drho + (1/(np.pi*e_base*delta_e))*(product_rule_rho)
        partials['CD','g'] = dCDdg =  dCD_0dg + (1/(np.pi*e_base*delta_e))*(product_rule_g)

        partials['CD','C_D0_base'] = dCDdCD0base =  dCD_0dCDbase 
        partials['CD','S_0'] = dCDdS0 =  dCD_0dS0 
        partials['CD','e_base'] = dCDdebase =  dCD_0debase - ((CL**2)/(np.pi*e_base*e_base*delta_e*inputs['AR']))
        partials['CD','delta_CD0'] = dCDddeltaCD0 =  dCD_0ddeltaCD0 
        partials['CD','delta_e'] = dCDddeltae =  dCD_0ddeltae - ((CL**2)/(np.pi*e_base*delta_e*delta_e*inputs['AR']))
        partials['CD','ks_base'] = dCDdks_base = delta_ks*(inputs['S']-S_0)
        partials['CD','delta_ks'] = dCDddelta_ks =  ks_base*(inputs['S']-S_0)

        partials['LD','V'] = (CD*dCLdV - CL*dCDdV)/(CD**2) 
        partials['LD','S'] = (CD*dCLdS - CL*dCDdS)/(CD**2)
        partials['LD','AR'] = (CD*dCLdAR - CL*dCDdAR)/(CD**2)
        partials['LD','m_total'] = (CD*dCLdm - CL*dCDdm)/(CD**2)
        partials['LD','rho'] = (CD*dCLdrho - CL*dCDdrho)/(CD**2)
        partials['LD','g'] = (CD*dCLdg - CL*dCDdg)/(CD**2)

        partials['LD','C_D0_base'] =  (0 - CL*dCDdCD0base)/(CD**2) 
        partials['LD','S_0'] = (0 - CL*dCDdS0)/(CD**2) 
        partials['LD','e_base'] = (0 - CL*dCDdebase)/(CD**2) 
        partials['LD','delta_CD0'] = (0 - CL*dCDddeltaCD0)/(CD**2) 
        partials['LD','delta_e'] = (0 - CL*dCDddeltae)/(CD**2) 
        partials['LD','ks_base'] = -(CL*dCDdks_base)/(CD**2)
        partials['LD','delta_ks'] = -(CL*dCDddelta_ks)/(CD**2)

        partials['WL','S'] = -(inputs['m_total']*g) / (inputs['S']**2)
        partials['WL','m_total'] = (g) / (inputs['S'])
        partials['WL','g'] = (inputs['m_total']) / (inputs['S'])


#function the returns a distribution of input variables on CI
def distribute_input(CI,base_val,sigma,n_points):
    mu = 1
    p_lower = (1-CI)/2 #lower end sample point
    p_upper = (1 + CI)/2

    p = np.linspace(p_lower,p_upper,n_points)

    delta_vec = erfinv(2*p - 1)*np.sqrt(2)*sigma + mu

    CDF_vec = (1.0/2.0)*erf((delta_vec-mu)/(np.sqrt(2)*sigma)) + 0.5

    u_vec = (delta_vec-mu)/(np.sqrt(2)*sigma)

    PDF_vec = (1.0/(sigma*np.sqrt(2*np.pi)))*np.exp(-(u_vec**2))

    #print(delta_vec)

    return delta_vec*base_val, CDF_vec, PDF_vec




def main():

    n_p = 5000

    prblm = om.Problem()
    prblm.model.add_subsystem('Aero',AeroDiscipline(vec_size=n_p))

    prblm.setup()

    prblm.set_val('Aero.g', 9.81)
    prblm.set_val('Aero.rho', 0.38)
    prblm.set_val('Aero.C_D0_base', parameters['CD0_base'])
    prblm.set_val('Aero.S_0', parameters['S_naught'])
    prblm.set_val('Aero.ks_base', parameters['ks_base'])
    prblm.set_val('Aero.e_base', parameters['e_oswald_base'])
    prblm.set_val('Aero.S', parameters['S_naught']*0.9)
    prblm.set_val('Aero.V', parameters['V_ref']) 
    prblm.set_val('Aero.AR', parameters['AR'])
    prblm.set_val('Aero.m_total', 80000.0)

    del_e, CDF_e, PDF_e = distribute_input(0.98,1.0,0.05,n_p)
    del_ks, CDF_ks, PDF_ks = distribute_input(0.98,1.0,0.15,n_p)
    del_CD0, CDF_CD0, PDF_CD0 = distribute_input(0.98,1.0,0.1,n_p)



    prblm.set_val('Aero.delta_CD0',del_CD0)
    prblm.set_val('Aero.delta_ks',del_ks)
    prblm.set_val('Aero.delta_e',del_e)

    

    prblm.run_model()

    CL = prblm.get_val('Aero.CL')

    CD = prblm.get_val('Aero.CD')

    LoD = prblm.get_val('Aero.LD')


    print("Expected Scalar CL:",CL,"\n")
    print("Expected Vector of L/D values\n",LoD)


    plt.rcParams.update({
        "text.usetex" : True,
        "font.family" : "serif"
    })

    figure_e, ax = plt.subplots(2,2)

    ax_pdf_de = ax[0,0]
    ax_cdf_de = ax[0,1]

    ax_pdf_e = ax[1,0]
    ax_cdf_e = ax[1,1]

    
    ax_pdf_de.plot(del_e,PDF_e,label=r"$\mathrm{PDF}(\delta_e)$")
    ax_pdf_de.legend()
    ax_pdf_de.set_xlabel(r"$\delta_e$")
    ax_pdf_de.set_ylabel("Probability Density")

    ax_cdf_de.plot(del_e,CDF_e,label=r"$\mathrm{CDF}(\delta_e)$")
    ax_cdf_de.legend()
    ax_cdf_de.set_xlabel(r"$\delta_e$")
    ax_cdf_de.set_ylabel("Cummulative Probability")

    ax_pdf_e.plot(del_e*parameters['e_oswald_base'],
                   PDF_e/(parameters['e_oswald_base']),
                   label=r"$\mathrm{PDF}(e_{\mathrm{oswald}})$")
    ax_pdf_e.legend()
    ax_pdf_e.set_xlabel(r"$e_{\mathrm{oswald}}$")
    ax_pdf_e.set_ylabel("Probability Density")

    ax_cdf_e.plot(del_e*parameters['e_oswald_base'],CDF_e,label=r"$\mathrm{CDF}(e_{\mathrm{oswald}})$")
    ax_cdf_e.legend()
    ax_cdf_e.set_xlabel(r"$e_{\mathrm{oswald}}$")
    ax_cdf_e.set_ylabel("Cummulative Probability")

    #ax_del_e.set_xlim([0,2])

    #plt.hist(del_e*0.8,bins=25)

    plt.show()

   




    




if __name__ == "__main__":
    main()