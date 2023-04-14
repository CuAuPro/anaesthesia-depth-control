function BIS = ContentrationToBIS(Uprop, Uremi, gamma, sigma, Emax)
   BIS = Emax - Emax*((Uprop)+(Uremi)+sigma.*(Uprop).*(Uremi)).^(gamma)./(1+(Uprop+Uremi+sigma.*Uprop.*Uremi).^(gamma));
end