function PAR_ppfd = PPFD_ChlF(wavelenth, radiance)

rad_umol=(3.14156.*radiance./100)./(10^(-6) * 6.022 * 10^(23) * 6.62 * 10^(-34) * 3 * 10^(8)./(wavelenth./10^9)) ;% (up_wl./10^9) from nm to meter       %
rad_umol_temp=nan(size(rad_umol));
for j=1:size(rad_umol,1)
    for ii=2:size(rad_umol,2)-1
        %  umol/m2/s/nm to  umol/m2/s
        rad_umol_temp(j,ii)=rad_umol(j,ii)*((wavelenth(1,ii+1)-wavelenth(1,ii-1))/2);% total chlf = radiance(x)*((wl(x+1)-wl(x-1))/2) x is wavelength unit mw/m2/sr
    end
end
% select the 400-700nm
PAR_wl=wavelenth>640 & wavelenth<850;
PAR_ppfd=sum(rad_umol_temp(:,PAR_wl),2,'omitnan'); %unit is  umol/m2/s
clear  HR_ref_PAR_wl

end