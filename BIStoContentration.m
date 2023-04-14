function Uprop = BIStoContentration(BIS, Uremi, gamma, sigma, Emax)
    Uprop = (((Emax-BIS)./BIS)^(1/gamma)-Uremi)/(1+sigma*Uremi);
end