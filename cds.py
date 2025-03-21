import bisect
import math
from enum import IntEnum

import numpy as np
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

# Product
class CDS:
    def __init__(self, maturity, freq, coupon, pay_accrued = True, r = None):
        self.maturity = maturity
        self.freq = freq
        self.coupon = coupon
        self.pay_accrued = pay_accrued
        self.r = r

# Curves
class Interpolator:
    def eval(self, t ):
        raise NotImplementedError

class PiecewiseFlatForwardInter(Interpolator):
    def __init__(self, dates, forwards):
        # assuming dates are sorted
        self.dates = dates
        self.forwards = forwards
        self.cums = [ 1 ]
        for n in range(1, len(dates)):
            self.cums.append(self.cums[n-1] *
                        math.exp(-self.forwards[n-1]*(self.dates[n] - self.dates[n-1])))

    def eval(self, t):
        left = bisect.bisect_left(self.dates, t)
        if left == len(self.dates):
            left -= 1

        ind = self.dates.index(left)
        return self.cums[ind] * math.exp(-self.forwards[ind]*(t-left))

class Curve:
    def __init__(self, inter):
        self.inter = inter

    def eval(self, t):
        return self.inter.eval(t)

class SurvCurve(Curve):
    def __init__(self, inter, r):
#        super(Curve, SurvCurve).__init__(inter)
        super().__init__(inter)
        self.r = r

# Pricers
class CdsPricingResults(IntEnum):
    NPV = 0
    LOSS = 1
    FIXED = 2
    ANNUITY = 3
    PAR_SPREAD = 4
    ACCRUED_INTEREST = 5 # not used

class CdsPricer:
    def __init__(self, disc_curve, surv_curve):
        self.disc_curve = disc_curve
        self.surve_curve = surv_curve

    def price(self, cds, buy_prot, notional):
        results = self.do_price_prot_buyer(cds) # assume returns array of results
        results = results * notional * (1 if buy_prot else -1)
        return results

    def do_price_prot_buyer(self, cds):
        raise NotImplementedError


class GenericCdsPricer(CdsPricer):
    def __init__(self, disc_curve, surv_curve, n_integ_step_cpn = 25):
#        super(CdsPricer, GenericCdsPricer).__init__(disc_curve, surv_curve)
        super().__init__(disc_curve, surv_curve)
        self.n_integ_step_cpn = n_integ_step_cpn

    def do_price_prot_buyer(self, cds):
        loss = 0
        coupon = 0
        acc_on_default = 0
        n_coupons = cds.maturity*cds.freq
        yf = 1/cds.freq
        step = yf/self.n_integ_step_cpn
        sp_prev = 1
        t = 0
        for n in range(n_coupons):
            coupon_start_t = t
            for integ_step in range(self.n_integ_step_cpn):
                t = t + step
                df = self.disc_curve.eval(t)
                sp = self.surve_curve.eval(t)
                ds = sp_prev - sp
                dfds = df * ds
                sp_prev = sp

                loss = loss + dfds
                acc_on_default = acc_on_default + (t-coupon_start_t) * dfds

            # the below could use the exact time,
            # but keeping it consistent with the accrued integration
            coupon = coupon + df*sp

        coupon = coupon * yf
        rec = self.surve_curve.r if cds.r is None else cds.r
        loss = loss * (1 - rec)
        annuity = coupon + (acc_on_default if cds.pay_accrued else 0)
        fixed = cds.coupon * annuity
        npv = loss - fixed
        par = loss/annuity
        accrued_interest = 0 # reserved for future use

        # Populate pricing results
        results = np.zeros(len(CdsPricingResults))
        results[CdsPricingResults.NPV] = npv
        results[CdsPricingResults.LOSS] =loss
        results[CdsPricingResults.FIXED] =fixed
        results[CdsPricingResults.PAR_SPREAD] =par
        results[CdsPricingResults.ANNUITY] =annuity
        results[CdsPricingResults.ACCRUED_INTEREST] =accrued_interest

        return results

# The below assumes that only a single degree of freedom is calibrated
# which is taken to be a single flat forward intensity

class FlatCalibration:
    def __init__(self, cds, r, disc_curve, npv = 0):
        self.cds = cds
        self.npv = npv
        self.inter = PiecewiseFlatForwardInter( [0], [0])
        self.surv_curve = SurvCurve(self.inter, r)
        self.pricer = GenericCdsPricer(disc_curve, self.surv_curve)

    def objective(self, x):
        self.inter.forwards[0] = x
        results = self.pricer.price(self.cds, True, 1)
        return results[CdsPricingResults.NPV] - self.npv

    def strip(self, bracket = [0, 100]):
        result = root_scalar( lambda x: self.objective(x), bracket = bracket )
        return self.surv_curve

# main code/tests

ir = 0.04
r = 0.4
c = 0.01
h = c / (1-r)
mat = 5
freq = 4

rates_interp = PiecewiseFlatForwardInter( [0], [ir])
disc_curve = Curve(rates_interp)

credit_interp = PiecewiseFlatForwardInter( [0], [h])
credit_curve = SurvCurve(credit_interp, r)

cds = CDS(mat, freq, c, True)

#pricing
print("Valuation")
pricer = GenericCdsPricer(disc_curve, credit_curve)
results = pricer.price(cds, True, 1.2)
for i in CdsPricingResults:
    print(i.name, " = ", results[i] )

#calibration
print("\nCalibration")
calib = FlatCalibration(cds, r, disc_curve, results[CdsPricingResults.NPV])
curve = calib.strip([0,2*h*(1-r)])
print("Hazard rate = ", curve.inter.forwards[0])

#recovery/npv graph
par =[]
rr = []

credit_interp = PiecewiseFlatForwardInter( [0], [h])
credit_curve = SurvCurve(credit_interp, r)
pricer = GenericCdsPricer(disc_curve, credit_curve)

for rec in range(0,10):
    cds.r = rec/10
    rr.append (cds.r)
    par.append(pricer.price(cds, True, 1)[CdsPricingResults.PAR_SPREAD])

plt.plot(rr, par)
plt.show()










