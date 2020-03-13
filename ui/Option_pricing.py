import sys
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5 import uic
import math

#import methods

sys.path.append("../model")
import black_scholes_cranknicolson as bscn
import black_scholes_implicit as bsi
import black_scholes_explicit as bse
import black_scholes_greeks as bsg


qtCreatorFile = "Option_pricing.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class Main(QMainWindow, Ui_MainWindow):          
    def __init__(self):
        super().__init__()
        self.setupUi(self)

        self.getModel.activated[str].connect(self.chooseMethod)


    def chooseMethod(self, text):
        cur_txt = text
        if cur_txt == 'Black Scholes - Crank Nicolson':
            self.calculate.clicked.connect(self.blackscholes_cranknicolson)
        elif cur_txt == 'Black Scholes - Implicit':
            self.calculate.clicked.connect(self.blackscholes_implicit)
        elif cur_txt == 'Black Scholes - Explicit':
            self.calculate.clicked.connect(self.blackscholes_explicit)
        elif cur_txt == 'Black Scholes - Greeks':
            self.calculate.clicked.connect(self.blackscholes_greeks)
        else:
            self.calculate.clicked.connect(self.blackscholes_cranknicolson)

  
    def blackscholes_implicit(self):
        #print("BS_Implicit")
        self.clearfield()
        S = float(self.stockPrice.text())
        K = float(self.exercisePrice.text())
        sigma = float(self.volatility.text())/100
        r = float(self.interestRate.text())/100
        q = float(self.yieldRate.text())/100
        dividendMethod = self.dividendMethod.text()
        valDate = self.valueDate.date().toPyDate()
        expDate = self.expirationDate.date().toPyDate()
        daysDiff = float((expDate - valDate).days)
        tau = daysDiff/365
        M = 100
        N = 100
     
        #Calling Black Scholes - Implicit Method -->
        bs_i_result = bsi.black_scholes_implicit(S,K,r,q,tau,sigma,M,N)
        
        # Printing onto UI
        theorecticalvalue_callString = str(round(bs_i_result[0],5)) #5 dp and convert to string
        self.theoreticalValue_Call.setText(theorecticalvalue_callString)

        theorecticalvalue_putString = str(round(bs_i_result[1],5)) #5 dp and convert to string
        self.theoreticalValue_Put.setText(theorecticalvalue_putString)

        


    def blackscholes_explicit(self):
        self.clearfield()
        S = float(self.stockPrice.text())
        K = float(self.exercisePrice.text())
        sigma = float(self.volatility.text())/100
        r = float(self.interestRate.text())/100
        q = float(self.yieldRate.text())/100
        dividendMethod = self.dividendMethod.text()
        valDate = self.valueDate.date().toPyDate()
        expDate = self.expirationDate.date().toPyDate()
        daysDiff = float((expDate - valDate).days)
        tau = daysDiff/365
        M = 100
        N = 1000
        
        #Calling Black Scholes - Explicit Method -->
        bs_e_result = bse.black_scholes_explicit(S,K,r,q,tau,sigma,M,N)
        
        # Printing onto UI
        theorecticalvalue_callString = str(round(bs_e_result[0],5)) #5 dp and convert to string
        self.theoreticalValue_Call.setText(theorecticalvalue_callString)

        theorecticalvalue_putString = str(round(bs_e_result[1],5)) #5 dp and convert to string
        self.theoreticalValue_Put.setText(theorecticalvalue_putString)


    def blackscholes_greeks(self):
        self.clearfield()
        S = float(self.stockPrice.text())
        K = float(self.exercisePrice.text())
        sigma = float(self.volatility.text())/100
        r = float(self.interestRate.text())/100
        q = float(self.yieldRate.text())/100
        dividendMethod = self.dividendMethod.text()
        valDate = self.valueDate.date().toPyDate()
        expDate = self.expirationDate.date().toPyDate()
        daysDiff = float((expDate - valDate).days)
        tau = daysDiff/365

        bs_g_tvalue = bsg.callPutOptionPrices(S,K,r,q,tau,sigma) #return [c,p]
        bs_g_deltas = bsg.deltas(S,K,r,q,tau,sigma) #return [deltaC, deltaP, deltaC100, deltaP100]
        bs_g_lambdas = bsg.lambdas(S,K,r,q,tau,sigma) #return [lambdasC, lambdasP]
        bs_g_bs = bsg.blackScholes(S,K,r,q,tau,sigma) #return (gammaC,gammaP,vegaC,vegaP,rhoC,rhoP,psiC,psiP,SSc,SSp,IVc,IVp,TVc,TVp,ZVc,ZVp)
        bs_g_gamma = bsg.gamma1Percent(S,K,r,q,tau,sigma) #return(GammaC,GammaP)
        bs_g_theta = bsg.theta1Day(S,K,r,q,tau,sigma) #return(thetaC1,thetaP1)
        bs_g_theta7 = bsg.theta7Days(S,K,r,q,tau,sigma) #return(thetaC7,thetaP7)

        #bs_g_tvalue
        theorecticalvalue_callString = str(round(bs_g_tvalue[0],5)) #5 dp and convert to string
        self.theoreticalValue_Call.setText(theorecticalvalue_callString)
        theorecticalvalue_putString = str(round(bs_g_tvalue[1],5)) #5 dp and convert to string
        self.theoreticalValue_Put.setText(theorecticalvalue_putString)

        #bs_g_deltas
        delta_callString = str(round(bs_g_deltas[0],5)) #5 dp and convert to string
        self.delta_call.setText(delta_callString)
        delta_putString = str(round(bs_g_deltas[1],5)) #5 dp and convert to string
        self.delta_put.setText(delta_putString)
        delta100_callString = str(round(bs_g_deltas[2],5)) #5 dp and convert to string
        self.delta100_call.setText(delta100_callString)
        delta100_putString = str(round(bs_g_deltas[3],5)) #5 dp and convert to string
        self.delta100_put.setText(delta100_putString)

        #bs_g_lambdas
        lambda_callString = str(round(bs_g_lambdas[0],5)) #5 dp and convert to string
        self.lambda_call.setText(lambda_callString)
        lambda_putString = str(round(bs_g_lambdas[1],5)) #5 dp and convert to string
        self.lambda_put.setText(lambda_putString)

        #bs_g_bs
        bsgbsList = ['gamma_call','gamma_put','vega_call','vega_put','rho_call','rho_put','psi_call','psi_put','strike_call','strike_put','intrinsic_call','intrinsic_put','timevalue_call','timevalue_put','zero_call','zero_put']

        for i in range(len(bsgbsList)):
            stringValue = str(round(bs_g_bs[i],5))
            x = getattr(self, bsgbsList[i])
            x.setText(stringValue)

        #bs_g_gamma
        gamma_callString = str(round(bs_g_gamma[0],5)) #5 dp and convert to string
        self.gamma1percent_call.setText(gamma_callString)
        gamma_putString = str(round(bs_g_gamma[1],5)) #5 dp and convert to string
        self.gamma1percent_put.setText(gamma_putString)

        #bs_g_theta
        theta_callString = str(round(bs_g_theta[0],5)) #5 dp and convert to string
        self.theta_call.setText(theta_callString)
        theta_putString = str(round(bs_g_theta[1],5)) #5 dp and convert to string
        self.theta_put.setText(theta_putString)

        #bs_g_theta7
        theta7_callString = str(round(bs_g_theta7[0],5)) #5 dp and convert to string
        self.theta7_call.setText(theta7_callString)
        theta7_putString = str(round(bs_g_theta7[1],5)) #5 dp and convert to string
        self.theta7_put.setText(theta7_putString)


    def blackscholes_cranknicolson(self):
        # print("BS_CN")
        self.clearfield()
        S = float(self.stockPrice.text())
        K = float(self.exercisePrice.text())
        sigma = float(self.volatility.text())/100
        r = float(self.interestRate.text())/100
        q = float(self.yieldRate.text())/100
        dividendMethod = self.dividendMethod.text()
        valDate = self.valueDate.date().toPyDate()
        expDate = self.expirationDate.date().toPyDate()
        daysDiff = float((expDate - valDate).days)
        tau = daysDiff/365
        M = 100
        N = 1000
       
        #Calling Black Scholes - Explicit Method -->
        bs_cn_result = bscn.blackScholes_cranknicolson(S,K,r,q,tau,sigma,M,N)
        
        # Printing onto UI
        theorecticalvalue_callString = str(round(bs_cn_result[0],5)) #5 dp and convert to string
        self.theoreticalValue_Call.setText(theorecticalvalue_callString)

        theorecticalvalue_putString = str(round(bs_cn_result[1],5)) #5 dp and convert to string
        self.theoreticalValue_Put.setText(theorecticalvalue_putString)


    def clearfield(self):
        self.theoreticalValue_Call.setText("")
        self.theoreticalValue_Put.setText("")
        self.delta_call.setText("")
        self.delta_put.setText("")
        self.delta100_call.setText("")
        self.delta100_put.setText("")
        self.lambda_call.setText("")
        self.lambda_put.setText("")
        self.gamma_call.setText("")
        self.gamma_put.setText("")
        self.gamma1percent_call.setText("")
        self.gamma1percent_put.setText("")
        self.theta_call.setText("")
        self.theta_put.setText("")
        self.theta7_call.setText("")
        self.theta7_put.setText("")
        self.strike_call.setText("")
        self.strike_put.setText("")
        self.intrinsic_call.setText("")
        self.intrinsic_put.setText("")
        self.timevalue_call.setText("")
        self.timevalue_put.setText("")
        self.zero_call.setText("")
        self.zero_put.setText("")      
        self.vega_call.setText("")
        self.vega_put.setText("")
        self.rho_call.setText("")
        self.rho_put.setText("")
        self.psi_call.setText("")
        self.psi_put.setText("")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main = Main()
    main.show()
    sys.exit(app.exec_())