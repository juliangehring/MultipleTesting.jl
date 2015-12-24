## lsl benchmark

using MultipleTesting

p = rand(1000);
m = 1000;

@time for i in 1:m
    pi0 = MultipleTesting.lsl_pi0(p)
end

@time for i in 1:m
    pi0 = MultipleTesting.lsl_pi0_vec(p)
end
