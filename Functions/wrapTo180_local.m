function ang = wrapTo180_local(ang)
ang = mod(ang + 180, 360) - 180;
end