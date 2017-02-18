module NK_baseline

    using NLsolve

    function steadystate(model) #,ys,exo)

        new_calibration = copy(model.calibration)

        for p in keys(model.calibration)
            v = model.calibration[p]
            expr = :($p=$v)
            eval(expr)
        end
        eta = model.calibration[:eta] # workaround

        PI=PIbar
        u=1
        q=1
        d=1
        phi=1
        m=0
        zeta=1
        mu_z=exp(LambdaYd)
        mu_I=exp(Lambdamu)
        mu_A=exp(LambdaA)

        # %set the parameter Lambdax
        Lambdax=mu_z

        # %set the parameter gammma1
        gammma1=mu_z*mu_I/betta-(1-delta)


        r=1*gammma1
        R=1+(PI*mu_z/betta-1)

        Rbar=R

        PIstar=((1-thetap*PI^(-(1-epsilon)*(1-chi)))/(1-thetap))^(1/(1-epsilon))
        PIstarw=((1-thetaw*PI^(-(1-chiw)*(1-eta))*mu_z^(-(1-eta)))/(1-thetaw))^(1/(1-eta))

        mc=(epsilon-1)/epsilon*(1-betta*thetap*PI^((1-chi)*epsilon))/(1-betta*thetap*PI^(-(1-epsilon)*(1-chi)))*PIstar
        w=(1-alppha)*(mc*(alppha/r)^alppha)^(1/(1-alppha))
        wstar=w*PIstarw
        vp=(1-thetap)/(1-thetap*PI^((1-chi)*epsilon))*PIstar^(-epsilon)
        vw=(1-thetaw)/(1-thetaw*PI^((1-chiw)*eta)*mu_z^eta)*PIstarw^(-eta)
        tempvaromega=alppha/(1-alppha)*w/r*mu_z*mu_I

        fobj(ld::Float64) = ( (1-betta*thetaw*mu_z^(eta-1)*PI^(-(1-chiw)*(1-eta)))/(1-betta*thetaw*mu_z^(eta*(1+gammma))*PI^(eta*(1-chiw)*(1+gammma)))
        -(eta-1)/eta*wstar/(varpsi*PIstarw^(-eta*gammma)*ld^gammma)*((1-h*mu_z^(-1))^(-1)-betta*h*(mu_z-h)^(-1))*
        ((mu_A*mu_z^(-1)*vp^(-1)*tempvaromega^alppha-tempvaromega*(1-(1-delta)*(mu_z*mu_I)^(-1)))*ld-vp^(-1)*Phi)^(-1) )
        fobj1(ld::Vector) = [fobj(ld[1])]

        sol = nlsolve(not_in_place(fobj1), [0.25])
        ld = sol.zero[1]

        l=vw*ld
        k=tempvaromega*ld
        x=(1-(1-delta)*(mu_z*mu_I)^(-1))*k
        yd=(mu_A/mu_z*k^alppha*ld^(1-alppha)-Phi)/vp
        c=(mu_A*mu_z^(-1)*vp^(-1)*tempvaromega^alppha-tempvaromega*(1-(1-delta)*(mu_z*mu_I)^(-1)))*ld-vp^(-1)*Phi
        lambda=(1-h*betta*mu_z^(-1))*(1-h/mu_z)^(-1)/c
        F=yd-1/(1-alppha)*w*ld
        f=(eta-1)/eta*wstar*PIstarw^(-eta)*lambda*ld/(1-betta*thetaw*mu_z^(eta-1)*PI^(-(1-chiw)*(1-eta)))
        f2=varpsi*d*phi*PIstarw^(-eta*(1+gammma))*ld^(1+gammma)/(1-betta*thetaw*(PI^chiw/PI)^(-eta*(1+gammma))*(wstar/wstar*mu_z)^(eta*(1+gammma)))

        g1=lambda*mc*yd/(1-betta*thetap*PI^((1-chi)*epsilon))
        g2=epsilon/(epsilon-1)*g1

        ######
        ######
        ######

        new_calibration[:PI] = PI
        new_calibration[:u] = u
        new_calibration[:q] = q
        new_calibration[:d] = d
        new_calibration[:phi] = phi
        new_calibration[:m] = m
        new_calibration[:zeta] = zeta
        new_calibration[:mu_z] = mu_z
        new_calibration[:mu_I] = mu_I
        new_calibration[:mu_A] = mu_A
        new_calibration[:Lambdax] = Lambdax
        new_calibration[:gammma1] = gammma1
        new_calibration[:r] = r
        new_calibration[:R] = R
        new_calibration[:Rbar] = Rbar
        new_calibration[:PIstar] = PIstar
        new_calibration[:PIstarw] = PIstarw
        new_calibration[:mc] = mc
        new_calibration[:w] = w
        new_calibration[:wstar] = wstar
        new_calibration[:vp] = vp
        new_calibration[:vw] = vw
        new_calibration[:tempvaromega] = tempvaromega
        new_calibration[:ld] = ld
        new_calibration[:l] = l
        new_calibration[:k] = k
        new_calibration[:x] = x
        new_calibration[:yd] = yd
        new_calibration[:c] = c
        new_calibration[:lambda] = lambda
        new_calibration[:F] = F
        new_calibration[:f] = f
        new_calibration[:f2] = f2
        new_calibration[:g1] = g1
        new_calibration[:g2] = g2

        return new_calibration

    end

    # not very satisfying and not to be imitated
    function update_calibration(model, new_calibration)
        for k in keys(new_calibration)
            model.calibration[k] = new_calibration[k]
        end
        for i=1:length(model.symbols[:endogenous])
            s = model.symbols[:endogenous][i]
            v = model.calibration[s]
            model.calibration_grouped[:endogenous][i] = v
        end
        for i=1:length(model.symbols[:parameters])
            s = model.symbols[:parameters][i]
            v = model.calibration[s]
            model.calibration_grouped[:parameters][i] = v
        end
    end

end
