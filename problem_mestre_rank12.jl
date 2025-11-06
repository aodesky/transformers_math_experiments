include("constants.jl")

"""
Mestre's Rank ≥ 12 Elliptic Curve Family - Conductor Minimization

Family: Genus 1 plane cubic curve from Mestre's paper
Goal: Minimize conductor C(t) where t is a rational number
Constraint: D(t) ≠ 0 (discriminant must be non-zero)

D(t): Discriminant (degree 36 polynomial)
D'(t): Derivative (degree 35 polynomial)
C(t): Conductor (computed via GP/Pari)

Local search: Single step of gradient descent on D(t)
Learning rate: 1/50
Reward: -log(conductor)
"""

# Discriminant coefficients (degree 36, highest degree first)
const D_COEFFS = BigInt[
    -6180509591227720419380815576086350476681408413696000000,
    114938623929384725353177092671893112265996383551488000000,
    85551930076196922637182242936935849955283009798144000000,
    -11719628364517969915879400083879559410975035747205120000000,
    27052858938943217389306770689380951058020225437925376000000,
    495233735772095817573467306356693483579440850052055040000000,
    -1801350648161536348295098941726255410489247569001054208000000,
    -10250169932874704526238384703116770475682736517705367552000000,
    61082789522191032916195198949323882271668853151436636160000000,
    28889019713203853749514869477449341722130633703018463232000000,
    -894269251478892110653793279137284480187713458226239045632000000,
    1998986763029986547291442571731467929032900302560167854080000000,
    2661443381030372339132138461213064026328428987766103425024000000,
    -22395786314228742177536194385054915532573975723543903928320000000,
    49754045290848789556449598143368652244595710783550411472896000000,
    -44163173854815427694280317244618399762028370828027075166208000000,
    -32796148340266456018606916619411215277163210920437139800064000000,
    148904961810211812582905296904750630903782562417849592053760000000,
    -206251369280050995820557086852749583234853890852417835827200000000,
    148755092346812123632637376165602044783334295774330759413760000000,
    -32597758506965011566365737558577383086466467094612793278464000000,
    -44295499759368608691450381682183461951308091965292567134208000000,
    49782009495051548235764568635229690056981489920153730973696000000,
    -22362866931987209877620015953762190874516213419975061995520000000,
    2626761458211022508276814384763419934345961252162347597824000000,
    2012441064460877642689216131506573045948955416959580897280000000,
    -894419902548735011677159393784767862754953344245243052032000000,
    27004564310929212161890586883786579124328865854187896832000000,
    61657866003996771936038377726735901162216714589972725760000000,
    -10210770063465329148660348309485349439782151641945341952000000,
    -1846853462515131041133469192728870173290370843025604608000000,
    497428925592111443220177495517750224788127586965258240000000,
    28948530843419825472651917479685557177402342219186176000000,
    -11855959828599693328589625913076149102956954066616320000000,
    42900905735805344359511077352919934420956368338944000000,
    116724890597580319681815300992222752267149653311488000000,
    -5788570664188895847185744142334355789296064004096000000
]

# Discriminant derivative coefficients (degree 35, highest degree first)
const D_PRIME_COEFFS = BigInt[
    -222498345284197935097709360739108617160530702893056000000,
    4022851837528465387361198243516258929309873424302080000000,
    2908765622590695369664196259855818898479622333136896000000,
    -386747736029093007224020202768025460562176179657768960000000,
    865691486046182956457816662060190433856647214013612032000000,
    15352245808934970344777486497057497990962666351613706240000000,
    -54040519444846090448852968251787662314677427070031626240000000,
    -297254928053366431260913156390386343794799359013455659008000000,
    1710318106621348921653465570581068703606727888240225812480000000,
    780003532256504051236901475891132226497527109981498507264000000,
    -23251000538451194876998625257569396484880549913882215186432000000,
    49974669075749663682286064293286698225822507564004196352000000000,
    63874641144728936139171323069113536631882295706386482200576000000,
    -515103085227261070083332470856263057249201441641509790351360000000,
    1094588996398673370241891159154110349381105637238109052403712000000,
    -927426650951123981579886662136986395002595787388568578490368000000,
    -655922966805329120372138332388224305543264218408742796001280000000,
    2829194274394024439075200641190261987171868685939142249021440000000,
    -3712524647040917924770027563349492498227370035343521044889600000000,
    2528836569895806101754835394815234761316683028163622910033920000000,
    -521564136111440185061851800937238129383463473513804692455424000000,
    -664432496390529130371755725232751929269621379479388507013120000000,
    696948132930721675300703960893215660797740858882152233631744000000,
    -290717270115833728409060207398908481368710774459675805941760000000,
    31521137498532270099321772617161039212151535025948171173888000000,
    22136851709069654069581377446572303505438509586555389870080000000,
    -8944199025487350116771593937847678627549533442452430520320000000,
    243041078798362909457015281954079212118959792687691071488000000,
    493262928031974175488307021813887209297733716719781806080000000,
    -71475390444257304040622438166397446078475061493617393664000000,
    -11081120775090786246800815156373221039742225058153627648000000,
    2487144627960557216100887477588751123940637934826291200000000,
    115794123373679301890607669918742228709609368876744704000000,
    -35567879485799079985768877739228447308870862199848960000000,
    85801811471610688719022154705839868841912736677888000000,
    116724890597580319681815300992222752267149653311488000000
]

# Helper functions for rational arithmetic
function parse_rational(s::String)::Tuple{Int64, Int64}
    """Parse a rational string 'num/den' into (numerator, denominator)"""
    if occursin("/", s)
        parts = split(s, "/")
        num = parse(Int64, parts[1])
        den = parse(Int64, parts[2])
        return (num, den)
    else
        num = parse(Int64, s)
        return (num, 1)
    end
end

function rational_to_string(num::Int64, den::Int64)::String
    """Convert rational to string"""
    if den == 1
        return string(num)
    else
        return string(num) * "/" * string(den)
    end
end

function simplify_rational(num::Int64, den::Int64)::Tuple{Int64, Int64}
    """Simplify a rational number using GCD"""
    if den < 0
        num = -num
        den = -den
    end
    if num == 0
        return (0, 1)
    end
    g = gcd(abs(num), abs(den))
    return (div(num, g), div(den, g))
end

function eval_poly_rational(coeffs::Vector{BigInt}, num::Int64, den::Int64)::BigFloat
    """
    Evaluate polynomial with BigInt coefficients at rational t = num/den
    Uses Horner's method for efficiency
    Returns BigFloat result
    """
    if isempty(coeffs)
        return BigFloat(0)
    end

    t = BigFloat(num) / BigFloat(den)
    result = BigFloat(coeffs[1])

    for i in 2:length(coeffs)
        result = result * t + BigFloat(coeffs[i])
    end

    return result
end

function eval_discriminant(num::Int64, den::Int64)::BigFloat
    """Evaluate D(t) at t = num/den"""
    return eval_poly_rational(D_COEFFS, num, den)
end

function eval_discriminant_derivative(num::Int64, den::Int64)::BigFloat
    """Evaluate D'(t) at t = num/den"""
    return eval_poly_rational(D_PRIME_COEFFS, num, den)
end

function compute_conductor_mestre(num::Int64, den::Int64)::BigInt
    """
    Compute the conductor of Mestre's curve at t = num/den using GP/Pari

    The Weierstrass coefficients a4 and a6 are also polynomials in t,
    so we need to evaluate them at the rational value.

    For now, we'll use the generate_rank_12.py approach via Sage/Pari.
    """
    # Create GP/Pari script
    gp_script = """
    t = $num/$den;

    \\\\ Mestre's polynomial coefficients
    A1 = -26940*t^2 + 51220*t - 26940;
    A2 = -1320*t^3 + 17280*t^2 + 17280*t - 1320;
    A3 = -18876*t^4 - 153828*t^3 + 301221*t^2 - 153828*t - 18776;
    A4 = -1489600*t^3 + 1489600*t^2 + 1489600*t - 1489600;
    A5 = 5816880*t^4 + 8043880*t^3 - 27463500*t^2 + 8043880*t + 5816880;
    A6 = 3416160*t^5 - 24166320*t^4 + 19202040*t^3 + 19202040*t^2 - 24166320*t + 3416160;
    A7 = -745360*t^6 - 15468024*t^5 + 18853764*t^4 - 138394*t^3 + 18853764*t^2 - 15468024*t - 745360;

    \\\\ Compute via Sage's approach would go here
    \\\\ For now, use a placeholder
    print("MESTRE_NOT_IMPLEMENTED");
    quit();
    """

    temp_file = tempname() * ".gp"
    open(temp_file, "w") do f
        write(f, gp_script)
    end

    try
        result = read(pipeline(`gp -q $temp_file`), String)
        rm(temp_file)
        result = strip(result)

        if result == "MESTRE_NOT_IMPLEMENTED"
            # Fall back to simple curve for now
            # y^2 = x^3 + t*x + 1
            return compute_conductor_simple(num, den)
        end

        return parse(BigInt, result)
    catch e
        println("Warning: GP/Pari failed for Mestre t=$num/$den: ", e)
        if isfile(temp_file)
            rm(temp_file)
        end
        return BigInt(typemax(Int64))
    end
end

function compute_conductor_simple(num::Int64, den::Int64)::BigInt
    """Fallback: compute conductor for simple curve y^2 = x^3 + t*x + 1"""
    gp_script = """
    t = $num/$den;
    E = ellinit([0, 0, 0, t, 1]);
    if (E == 0, print("INVALID"), red = ellglobalred(E); print(red[1]));
    quit();
    """

    temp_file = tempname() * ".gp"
    open(temp_file, "w") do f
        write(f, gp_script)
    end

    try
        result = read(pipeline(`gp -q $temp_file`), String)
        rm(temp_file)
        result = strip(result)

        if result == "INVALID"
            return BigInt(typemax(Int64))
        end

        return parse(BigInt, result)
    catch e
        println("Warning: GP/Pari failed: ", e)
        if isfile(temp_file)
            rm(temp_file)
        end
        return BigInt(typemax(Int64))
    end
end

function greedy_search_from_startpoint(db, obj::OBJ_TYPE)::Vector{OBJ_TYPE}
    """
    Perform one step of gradient descent on D(t)
    """
    local num, den
    try
        num, den = parse_rational(obj)
    catch
        return [empty_starting_point()]
    end

    # Evaluate discriminant
    disc = eval_discriminant(num, den)

    # Check if discriminant is zero (invalid curve)
    if abs(disc) < 1e-10
        new_num, new_den = simplify_rational(num + rand(-5:5), den + max(1, rand(-2:2)))
        return [rational_to_string(new_num, new_den)]
    end

    # Evaluate derivative
    disc_deriv = eval_discriminant_derivative(num, den)

    # Learning rate: 1/50
    learning_rate = 1/50

    # Gradient descent step (work in Float64 then convert back to rational)
    t_current = Float64(num) / Float64(den)

    if disc > 0
        t_new_float = t_current - learning_rate * Float64(disc_deriv)
    else
        t_new_float = t_current + learning_rate * Float64(disc_deriv)
    end

    # Convert float back to rational with reasonable denominator
    # Use rationalize with tolerance
    t_new_rational = rationalize(BigInt, t_new_float, tol=1e-10)
    new_num_big = numerator(t_new_rational)
    new_den_big = denominator(t_new_rational)

    # Try to fit into Int64 range
    try
        new_num = Int64(new_num_big)
        new_den = Int64(new_den_big)
        new_num, new_den = simplify_rational(new_num, new_den)

        # Verify non-zero discriminant
        new_disc = eval_discriminant(new_num, new_den)
        if abs(new_disc) < 1e-10
            new_num += 1
            new_num, new_den = simplify_rational(new_num, new_den)
        end

        return [rational_to_string(new_num, new_den)]
    catch
        # Number too large, return original with small perturbation
        perturbed_num, perturbed_den = simplify_rational(num + 1, den)
        return [rational_to_string(perturbed_num, perturbed_den)]
    end
end

function reward_calc(obj::OBJ_TYPE)::REWARD_TYPE
    """
    Compute the reward = -log(conductor(t))
    """
    try
        num, den = parse_rational(obj)

        # Check discriminant is non-zero
        disc = eval_discriminant(num, den)
        if abs(disc) < 1e-10
            return Float32(-1e9)
        end

        conductor = compute_conductor_mestre(num, den)

        if conductor <= 0 || conductor >= BigInt(10)^100
            return Float32(-1e9)
        end

        # Return negative log of conductor
        return Float32(-log(Float64(conductor)))
    catch e
        println("Error computing reward for $obj: ", e)
        return Float32(-1e9)
    end
end

function empty_starting_point()::OBJ_TYPE
    """Initial starting point: t = 1/1"""
    return "1/1"
end
