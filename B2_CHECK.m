function val = B2_CHECK(F2,v2)
    if v2 <= 0
        val = 0.1;
    else
        B_guess = 0.6 * F2;
        val = B_guess / v2;
    end
end
