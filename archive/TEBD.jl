using ITensors, ITensorMPS

t_hop  = 1.0
V      = 1.0
tp_hop = 0.1
Vp     = 0.1
μ      = 0.0

let
    N      = 100
    cutoff = 1E-8
    maxdim = 32
    tau    = 0.1
    ttotal = 5.0

    s = siteinds("S=1/2", N; conserve_qns=true)

    h_eff = -μ + V + Vp  # effective on-site field

    function make_NN_gate(j, tau)
        s1, s2 = s[j], s[j+1]
        hj = -t_hop   * op("S+",s1) * op("S-",s2) +
             -t_hop   * op("S-",s1) * op("S+",s2) +
              V       * op("Sz",s1) * op("Sz",s2) +
             -(h_eff/2) * op("Sz",s1) * op("Id",s2) +
             -(h_eff/2) * op("Id",s1) * op("Sz",s2)
        return exp(-im * tau * hj)
    end

    function make_NNN_gate(j, tau)
        s1, s2 = s[j], s[j+2]
        hj = -tp_hop * op("S+",s1) * op("S-",s2) +
             -tp_hop * op("S-",s1) * op("S+",s2) +
              Vp     * op("Sz",s1) * op("Sz",s2)
        return exp(-im * tau * hj)
    end

    gates_NN_odd  = [make_NN_gate(j, tau) for j in 1:2:(N-1)]
    gates_NN_even = [make_NN_gate(j, tau) for j in 2:2:(N-1)]
    gates_NNN_L1  = [make_NNN_gate(j, tau) for j in 1:(N-2) if (j-1) % 4 ∈ (0,1)]
    gates_NNN_L2  = [make_NNN_gate(j, tau) for j in 1:(N-2) if (j-1) % 4 ∈ (2,3)]

    gates = vcat(gates_NN_odd, gates_NN_even, gates_NNN_L1, gates_NNN_L2)

    psi = MPS(s, n -> isodd(n) ? "Up" : "Dn")
    c   = div(N, 2)

    for t in 0.0:tau:ttotal
        Sz = expect(psi, "Sz"; sites=c)
        nj = Sz + 0.5
        println("t = $t,  <n_$c> = $nj")

        t ≈ ttotal && break

        psi = apply(gates, psi; cutoff, maxdim)
        normalize!(psi)
    end
    return
end