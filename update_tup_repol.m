function [tup, trepol, beat_num] = update_tup_repol(ti, dt, Vm, Vm_old, Vthresh, ...
    tup, trepol, beat_num)
iact = find((Vm>Vthresh).*(Vm_old<Vthresh));
for j = 1:length(iact)
    y1 = Vm_old(iact(j)); y2 = Vm(iact(j));
    m = (y2-y1)/dt;
    tup(iact(j),beat_num(iact(j))) = ti - (y1-Vthresh)./m;
    %         save(times_name,'tup','trepol');
end
irepol = find((Vm<Vthresh).*(Vm_old>Vthresh));
for j = 1:length(irepol)
    if ti-tup(irepol(j),beat_num(irepol(j)))>.1
        y1 = Vm_old(irepol(j)); y2 = Vm(irepol(j));
        m = (y2-y1)/dt;
        trepol(irepol(j),beat_num(irepol(j))) = ti - (y1-Vthresh)./m;
        beat_num(irepol(j)) = beat_num(irepol(j)) + 1;
    end
    %         save(times_name,'tup','trepol');
end

end

