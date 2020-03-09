import sys
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5 import uic
import math

#import methods
# import binomial as bi
# import montecarlo as mc
# import blackscholes as bs

qtCreatorFile = "Option_pricing.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class Main(QMainWindow, Ui_MainWindow):          
    def __init__(self):
        super().__init__()
        self.setupUi(self)

        self.model.activated[str].connect(self.chooseMethod)


    def chooseMethod(self, text):
        cur_txt = text
        # if  cur_txt == 'Binomial':
        #     self.calculate.clicked.connect(self.binomial)
        # elif cur_txt == 'Monte Carlo Simulation':
        #     self.calculate.clicked.connect(self.montecarlo)
        # elif cur_txt == 'Black Scholes':
        #     self.calculate.clicked.connect(self.blackscholes)
        # elif cur_txt == 'Black Scholes - Implicit':
        #     self.calculate.clicked.connect(self.blackscholes_implicit)
        # elif cur_txt == 'Black Scholes - Explicit':
        #     self.calculate.clicked.connect(self.blackscholes_explicit)
        # elif cur_txt == 'Black Scholes - Crank Nicolson':
        #     self.calculate.clicked.connect(self.blackscholes_cranknicolson)
        # else:
        #     self.calculate.clicked.connect(self.blackscholes)

  
        
    def PB_C(self):
        #--------------------START------------------------------
        print('Write code here.')
        print(self.lineEdit_LoanAmount.text(), 
              type(self.lineEdit_LoanAmount.text()))
        print(self.lineEdit_RepaymentPeriod.text())
        print(self.lineEdit_InterestRateOfLoan.text())
        self.lineEdit_MonthlyInstallment.setText('Write Answer')
        #-------------------- END ------------------------------
        # 0. Comment the above 6 lines of code.
        # 1. Use the template to get strings from GUI components.        
        # 2. Convert strings to numbers.
        # 3. Use variables.
        # 4. Write an expression to calcualte, including rounding.
        # 5. Convert the result from number to string.
        # 6. Use the template to output the string to the GUI component.      

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main = Main()
    main.show()
    sys.exit(app.exec_())