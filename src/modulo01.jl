# Titulos constituintes dos indices S&P500, DJIA, NIKKEY225, FTSE100, DAX30, CAC40, PSI20 e EURO STOXX 50, em 31/03/2021

# The S&P 500 stock market index, maintained by S&P Dow Jones Indices, comprises 505 common stocks issued by 500 large-cap companies and traded on American stock exchanges (including the 30 companies that compose the Dow Jones Industrial Average), and covers about 80 percent of the American equity market by capitalization. Although called the S&P 500, the index contains 505 stocks because it includes two share classes of stock from 5 of its component companies.

const SP500n = ["A", "AAL", "AAP", "AAPL", "ABBV", "ABC", "ABMD", "ABT", "ACN", "ADBE", "ADI", "ADM", "ADP", "ADSK", "AEE", "AEP", "AES", "AFL", "AIG", "AIZ", "AJG", "AKAM", "ALB", "ALGN", "ALK", "ALL", "ALLE", "ALXN", "AMAT", "AMCR", "AMD", "AME", "AMGN", "AMP", "AMT", "AMZN", "ANET", "ANSS", "ANTM", "AON", "AOS", "APA", "APD", "APH", "APTV", "ARE", "ATO", "ATVI", "AVB", "AVGO", "AVY", "AWK", "AXP", "AZO", "BA", "BAC", "BAX", "BBY", "BDX", "BEN", "BF-B", "BIIB", "BIO", "BK", "BKNG", "BKR", "BLK", "BLL", "BMY", "BR", "BRK-B", "BSX", "BWA", "BXP", "C", "CAG", "CAH", "CARR", "CAT", "CB", "CBOE", "CBRE", "CCI", "CCL", "CDNS", "CDW", "CE", "CERN", "CF", "CFG", "CHD", "CHRW", "CHTR", "CI", "CINF", "CL", "CLX", "CMA", "CMCSA", "CME", "CMG", "CMI", "CMS", "CNC", "CNP", "COF", "COG", "COO", "COP", "COST", "CPB", "CPRT", "CRM", "CSCO", "CSX", "CTAS", "CTLT", "CTSH", "CTVA", "CTXS", "CVS", "CVX", "CZR", "D", "DAL", "DD", "DE", "DFS", "DG", "DGX", "DHI", "DHR", "DIS", "DISCA", "DISCK", "DISH", "DLR", "DLTR", "DOV", "DOW", "DPZ", "DRE", "DRI", "DTE", "DUK", "DVA", "DVN", "DXC", "DXCM", "EA", "EBAY", "ECL", "ED", "EFX", "EIX", "EL", "EMN", "EMR", "ENPH", "EOG", "EQIX", "EQR", "ES", "ESS", "ETN", "ETR", "ETSY", "EVRG", "EW", "EXC", "EXPD", "EXPE", "EXR", "F", "FANG", "FAST", "FB", "FBHS", "FCX", "FDX", "FE", "FFIV", "FIS", "FISV", "FITB", "FLIR", "FLT", "FMC", "FOX", "FOXA", "FRC", "FRT", "FTNT", "FTV", "GD", "GE", "GILD", "GIS", "GL", "GLW", "GM", "GNRC", "GOOG", "GOOGL", "GPC", "GPN", "GPS", "GRMN", "GS", "GWW", "HAL", "HAS", "HBAN", "HBI", "HCA", "HD", "HES", "HFC", "HIG", "HII", "HLT", "HOLX", "HON", "HPE", "HPQ", "HRL", "HSIC", "HST", "HSY", "HUM", "HWM", "IBM", "ICE", "IDXX", "IEX", "IFF", "ILMN", "INCY", "INFO", "INTC", "INTU", "IP", "IPG", "IPGP", "IQV", "IR", "IRM", "ISRG", "IT", "ITW", "IVZ", "J", "JBHT", "JCI", "JKHY", "JNJ", "JNPR", "JPM", "K", "KEY", "KEYS", "KHC", "KIM", "KLAC", "KMB", "KMI", "KMX", "KO", "KR", "KSU", "L", "LB", "LDOS", "LEG", "LEN", "LH", "LHX", "LIN", "LKQ", "LLY", "LMT", "LNC", "LNT", "LOW", "LRCX", "LUMN", "LUV", "LVS", "LW", "LYB", "LYV", "MA", "MAA", "MAR", "MAS", "MCD", "MCHP", "MCK", "MCO", "MDLZ", "MDT", "MET", "MGM", "MHK", "MKC", "MKTX", "MLM", "MMC", "MMM", "MNST", "MO", "MOS", "MPC", "MPWR", "MRK", "MRO", "MS", "MSCI", "MSFT", "MSI", "MTB", "MTD", "MU", "MXIM", "NCLH", "NDAQ", "NEE", "NEM", "NFLX", "NI", "NKE", "NLOK", "NLSN", "NOC", "NOV", "NOW", "NRG", "NSC", "NTAP", "NTRS", "NUE", "NVDA", "NVR", "NWL", "NWS", "NWSA", "NXPI", "O", "ODFL", "OKE", "OMC", "ORCL", "ORLY", "OTIS", "OXY", "PAYC", "PAYX", "PBCT", "PCAR", "PEAK", "PEG", "PENN", "PEP", "PFE", "PFG", "PG", "PGR", "PH", "PHM", "PKG", "PKI", "PLD", "PM", "PNC", "PNR", "PNW", "POOL", "PPG", "PPL", "PRGO", "PRU", "PSA", "PSX", "PVH", "PWR", "PXD", "PYPL", "QCOM", "QRVO", "RCL", "RE", "REG", "REGN", "RF", "RHI", "RJF", "RL", "RMD", "ROK", "ROL", "ROP", "ROST", "RSG", "RTX", "SBAC", "SBUX", "SCHW", "SEE", "SHW", "SIVB", "SJM", "SLB", "SNA", "SNPS", "SO", "SPG", "SPGI", "SRE", "STE", "STT", "STX", "STZ", "SWK", "SWKS", "SYF", "SYK", "SYY", "T", "TAP", "TDG", "TDY", "TEL", "TER", "TFC", "TFX", "TGT", "TJX", "TMO", "TMUS", "TPR", "TRMB", "TROW", "TRV", "TSCO", "TSLA", "TSN", "TT", "TTWO", "TWTR", "TXN", "TXT", "TYL", "UA", "UAA", "UAL", "UDR", "UHS", "ULTA", "UNH", "UNM", "UNP", "UPS", "URI", "USB", "V", "VAR", "VFC", "VIAC", "VLO", "VMC", "VNO", "VRSK", "VRSN", "VRTX", "VTR", "VTRS", "VZ", "WAB", "WAT", "WBA", "WDC", "WEC", "WELL", "WFC", "WHR", "WLTW", "WM", "WMB", "WMT", "WRB", "WRK", "WST", "WU", "WY", "WYNN", "XEL", "XLNX", "XOM", "XRAY", "XYL", "YUM", "ZBH", "ZBRA", "ZION", "ZTS", "^GSPC"]


const DJ30n = [:AAPL, :AMGN, :AXP, :BA, :CAT, :CRM, :CSCO, :CVX, :DD, :DIS, :GS, :HD, :HON, :IBM, :INTC, :JNJ, :JPM, :KO, :MCD, :MMM, :MRK,:MSFT, :NKE, :PG, :TRV, :UNH, :V,  :VZ, :WBA, :WMT, "^DJI", "^GSPC"]


const NIKKEI225n = ["6857.T", "8267.T", "5201.T", "2802.T", "6770.T", "6113.T", "9202.T", "8304.T", "2502.T", "3407.T", "4503.T", "7832.T", "5108.T", "7751.T", "6952.T", "9022.T", "9502.T", "4519.T", "7762.T", "1721.T", "7186.T", "8253.T", "4751.T", "7912.T", "8750.T", "4568.T", "6367.T", "1925.T", "8601.T", "2432.T", "4061.T", "6902.T", "4324.T", "4631.T", "5714.T", "9020.T", "6361.T", "4523.T", "5020.T", "6954.T", "9983.T", "6504.T", "4901.T", "5803.T", "6702.T", "8354.T", "5801.T", "6674.T", "1808.T", "7205.T", "6305.T", "7004.T", "6501.T", "7267.T", "5019.T", "7013.T", "1605.T", "3099.T", "7202.T", "8001.T", "3086.T", "8697.T", "6178.T", "2914.T", "5411.T", "1963.T", "6473.T", "1812.T", "4452.T", "7012.T", "9107.T", "9433.T", "9008.T", "9009.T", "2801.T", "2503.T", "5406.T", "6301.T", "9766.T", "4902.T", "6326.T", "3405.T", "6971.T", "4151.T", "2413.T", "8002.T", "1333.T", "8252.T", "8628.T", "7261.T", "2269.T", "6479.T", "4188.T", "8058.T", "6503.T", "8802.T", "7011.T", "9301.T", "5711.T", "7211.T", "8306.T", "8031.T", "4183.T", "7003.T", "8801.T", "5706.T", "9104.T", "8411.T", "8725.T", "6701.T", "3659.T", "5333.T", "2282.T", "2871.T", "7731.T", "5214.T", "9062.T", "5703.T", "3863.T", "5202.T", "5401.T", "1332.T", "9432.T", "9101.T", "4021.T", "7201.T", "2002.T", "3105.T", "6988.T", "8604.T", "6471.T", "6472.T", "9613.T", "1802.T", "9007.T", "3861.T", "6703.T", "6103.T", "7733.T", "6645.T", "9532.T", "4578.T", "5541.T", "6752.T", "4755.T", "6098.T", "8308.T", "7752.T", "2501.T", "7735.T", "9735.T", "6724.T", "1928.T", "3382.T", "6753.T", "1803.T", "4063.T", "8303.T", "4507.T", "4911.T", "4004.T", "9412.T", "9434.T", "9984.T", "2768.T", "8630.T", "6758.T", "7270.T", "3436.T", "4005.T", "8053.T", "4506.T", "5802.T", "6302.T", "5713.T", "8316.T", "8309.T", "5232.T", "8830.T", "7269.T", "8795.T", "5233.T", "1801.T", "6976.T", "2531.T", "8233.T", "4502.T", "6762.T", "3401.T", "4543.T", "8331.T", "5631.T", "9503.T", "8355.T", "5101.T", "9001.T", "9602.T", "5707.T", "5301.T", "8766.T", "4043.T", "9501.T", "8035.T", "9531.T", "8804.T", "9005.T", "3289.T", "7911.T", "3402.T", "4042.T", "5332.T", "5901.T", "3101.T", "7203.T", "8015.T", "4704.T", "4208.T", "3103.T", "9021.T", "7951.T", "7272.T", "9064.T", "6506.T", "6841.T", "4689.T", "^N225"]


const FTSE100n = ["III.L", "ADM.L", "AAL.L", "ANTO.L", "AHT.L", "ABF.L", "AZN.L", "AUTO.L", "AVST.L", "AVV.L", "AV.L", "BME.L", "BA.L", "BARC.L", "BDEV.L", "BKG.L", "BHP.L", "BP.L", "BATS.L", "BLND.L", "BT-A.L", "BNZL.L", "BRBY.L", "CCH.L", "CPG.L", "CRH.L", "CRDA.L", "DCC.L", "DGE.L", "ENT.L", "EVR.L", "EXPN.L", "FERG.L", "FLTR.L", "FRES.L", "GSK.L", "GLEN.L", "HLMA.L", "HL.L", "HIK.L", "HSBA.L", "IHG.L", "IMB.L", "INF.L", "ICP.L", "IAG.L", "ITRK.L", "JD.L", "JMAT.L", "JET.L", "KGF.L", "LAND.L", "LGEN.L", "LLOY.L", "LSEG.L", "MNG.L", "MRO.L", "MNDI.L", "NG.L", "NWG.L", "NXT.L", "OCDO.L", "PSON.L", "PSH.L", "PSN.L", "PHNX.L", "POLY.L", "PRU.L", "RB.L", "REL.L", "RSW.L", "RTO.L", "RMV.L", "RIO.L", "RR.L", "RDSA.L", "RSA.L", "SGE.L", "SBRY.L", "SDR.L", "SMT.L", "SGRO.L", "SVT.L", "SMDS.L", "SMIN.L", "SN.L", "SKG.L", "SPX.L", "SSE.L", "STAN.L", "SLA.L", "STJ.L", "TW.L", "TSCO.L", "ULVR.L", "UU.L", "VOD.L", "WEIR.L", "WTB.L", "WPP.L", "^FTSE"]


const DAX30n = ["ADS.DE", "ALV.DE", "BAS.DE", "BAYN.DE", "BEI.DE", "BMW.DE", "CON.DE", "1COV.DE", "DAI.DE", "DHER.DE", "DBK.DE", "DB1.DE", "DPW.DE", "DTE.DE", "DWNI.DE", "EOAN.DE", "FRE.DE", "FME.DE", "HEI.DE", "HEN3.DE", "IFX.DE", "LIN.DE", "MRK.DE", "MTX.DE", "MUV2.DE", "RWE.DE", "SAP.DE", "SIE.DE", "VOW3.DE", "VNA.DE", "^GDAXI"]


const CAC40n = ["AI.PA", "AIR.PA", "ALO.PA", "MT.AS", "ATO.PA", "CS.PA", "BNP.PA", "EN.PA", "CAP.PA", "CA.PA", "ACA.PA", "BN.PA", "DSY.PA", "ENGI.PA", "EL.PA", "RMS.PA", "KER.PA", "OR.PA", "LR.PA", "MC.PA", "ML.PA", "ORA.PA", "RI.PA", "PUB.PA", "RNO.PA", "SAF.PA", "SGO.PA", "SAN.PA", "SU.PA", "GLE.PA", "STLA.PA", "STM.PA", "TEP.PA", "HO.PA", "FP.PA", "URW.AS", "VIE.PA", "DG.PA", "VIV.PA", "WLN.PA", "^FCHI"]


const FTSEMIBn = ["A2A.MI", "AMP.MI", "ATL.MI", "AZM.MI", "BGN.MI", "BAMI.MI", "BPE.MI", "BRE.MI", "BZU.MI", "CPR.MI", "CNHI.MI", "DIA.MI", "ENEL.MI", "ENI.MI", "EXO.MI", "RACE.MI", "FBK.MI", "G.MI", "ISP.MI", "IG.MI", "JUVE.MI", "LDO.MI", "MB.MI", "MONC.MI", "PIRC.MI", "PST.MI", "PRY.MI", "REC.MI", "SPM.MI", "SFER.MI", "SRG.MI", "STLA.MI", "STM.MI", "TIT.MI", "TEN.MI", "TRN.MI", "UBI.MI", "UCG.MI", "UNI.MI", "US.MI", "FTSEMIB.MI"]


const PSI20n = ["ALTR.LS", "BCP.LS", "COR.LS", "CTT.LS", "EDP.LS", "EDPR.LS", "GALP.LS", "IBS.LS", "JMT.LS", "EGL.LS", "NBA.LS", "NOS.LS", "NVG.LS", "PHR.LS", "RAM.LS", "RENE.LS", "SEM.LS", "SON.LS", "PSI20.LS"]


const EUROSTOXX50n = ["ADS.DE", "AD.AS", "AI.PA", "AIR.PA", "ALV.DE", "ABI.BR", "ASML.AS", "AMS.MC", "CS.PA", "BBVA.MC", "SAN.MC", "BAS.DE", "BAYN.DE", "BMW.DE", "BNP.PA", "CRG.IR", "SGO.PA", "DAI.DE", "DPW.DE", "DTE.DE", "ENEL.MI", "ENGI.PA", "ENI.MI", "EOAN.DE", "EL.PA", "FRE.DE", "BN.PA", "IBE.MC", "ITX.MC", "INGA.AS", "ISP.MI", "LIN.DE", "OR.PA", "MC.PA", "MUV2.DE", "NOKIA.HE", "ORA.PA", "PHIA.AS", "SAF.PA", "SAN.PA", "SAP.DE", "SU.PA", "SIE.DE", "GLE.PA", "TEF.MC", "FP.PA", "URW.AS", "UNA.AS", "DG.PA", "VIV.PA", "VOW.DE"]


const ASX50n = ["A2M.AX", "AGL.AX", "AMC.AX", "ANZ.AX", "APA.AX", "ALD.AX", "ALL.AX", "ASX.AX", "AZJ.AX", "BHP.AX", "BXB.AX", "COH.AX", "COL.AX", "CBA.AX", "CPU.AX", "CSL.AX", "DXS.AX", "FMG.AX", "GMG.AX", "GPT.AX", "IAG.AX", "JHX.AX", "LLC.AX", "MQG.AX", "MPL.AX", "MGR.AX", "NAB.AX", "NCM.AX", "OSH.AX", "ORI.AX", "ORG.AX", "QAN.AX", "QBE.AX", "RHC.AX", "RIO.AX", "STO.AX", "SCG.AX", "SHL.AX", "S32.AX", "SGP.AX", "SUN.AX", "SYD.AX", "TLS.AX", "TCL.AX", "TWE.AX", "VCX.AX", "WES.AX", "WBC.AX", "WPL.AX", "WOW.AX", "^AFLI"] 


const HKSEn = ["0005.HK", "0011.HK", "0388.HK", "0939.HK", "1299.HK", "1398.HK", "2318.HK", "2388.HK", "2628.HK", "3328.HK", "3988.HK", "0002.HK", "0003.HK", "0006.HK", "1038.HK", "0012.HK", "0016.HK", "0017.HK", "0101.HK", "0688.HK", "0823.HK", "0960.HK", "1109.HK", "1113.HK", "1997.HK", "2007.HK", "0001.HK", "0027.HK", "0066.HK", "0175.HK", "0241.HK", "0267.HK", "0288.HK", "0386.HK", "0669.HK", "0700.HK", "0762.HK", "0857.HK", "0883.HK", "0941.HK", "1044.HK", "1093.HK", "1177.HK", "1810.HK", "1876.HK", "1928.HK", "2018.HK", "2020.HK", "2269.HK", "2313.HK", "2319.HK", "2382.HK", "3690.HK", "6862.HK", "9988.HK", "^HSI"]

const SMIn = ["NESN.SW", "NOVN.SW", "ROG.SW", "ZURN.SW", "UBSG.SW", "ABBN.SW", "CFR.SW", "LONN.SW", "SIKA.SW", "ALC.SW", "GIVN.SW", "LHN.SW", "SCMN.SW", "CSGN.SW", "PGHN.SW", "SREN.SW", "GEBN.SW", "SGSN.SW", "SLHN.SW", "UHR.SW"]

const FXn = ["EURUSD=X", "EURGBP=X", "EURJPY=X", "EURCHF=X", "EURSEK=X", "EURDKK=X", "EURCAD=X", "EURAUD=X", "EURKRW=X", "HKD=X"]



"""
    yahoo(simbolos, data_inicial, data_final, intervalo = "1d")
    
Obter dados de Yahoo.Finance

## Argumentos
* `simbolos`: simbolos de mercado, e.g "AAPL", "^GSPC", ["IBM", "FB"]
* `data_inicial`: Data inicial do periodo, e.g. Date(2019,12,20) ou DateTime(2019,12,20,8,30,0)
* `data_final`: Data final do periodo, e.g. Date(2020,12,20) ou DateTime(2020,12,20,8,30,0)
* `intervalo`: Intervalo de amostragem, e.g. "1d" (default), "1wk", "1mo", "3mo"
"""
function yahoo(symbol, date1 = Date(1900,1,1), date2 = Date(Dates.now()), interval::String = "1d")
    if isa(symbol, String) == true
        nsymb = 1
    elseif isa(symbol, Symbol) == true 
        nsymb = 1
    else
        nsymb = length(symbol)
    end
    
    date1 = DateTime(date1)
    date2 = DateTime(date2)
    if date1 > date2
        date1, date2 = date2, date1
    end
    date2 = date2 +  Day(1)  
    from = string(round(Int64, datetime2unix(date1)))
    to = string(round(Int64, datetime2unix(date2)))
    host = rand(["query1", "query2"])
    print("Concluído:")
    for i in 1:nsymb
        if nsymb == 1 && (isa(symbol, String) == true || isa(symbol, Symbol) == true) 
            symb = String(symbol)
        else
            symb = String(symbol[i])
        end
        print(" ", symb)
        url = "https://$host.finance.yahoo.com/v7/finance/chart/$symb?&interval=$interval&period1=$from&period2=$to"
        response = HTTP.get(url, cookies = true)
        body = JSON.parse(String(response.body))["chart"]["result"][1]   
        # values = body["indicators"]["quote"][1]
        a = replace!(body["indicators"]["adjclose"][1]["adjclose"], nothing => missing)
        new = DataFrame(
            Time = Dates.Date.(unix2datetime.(body["timestamp"])),
            # Open = convert(Vector{Float64},values["open"]),
            # High =  convert(Vector{Float64},values["high"]),
            # Low = convert(Vector{Float64},values["low"]),
            # Close = convert(Vector{Float64},values["close"]),
            #AdjClose = convert(Vector{Float64},body["indicators"]["adjclose"][1]["adjclose"]),
            AdjClose = convert(Vector{Union{Missing, Float64}},a),
            # Volume = convert(Vector{Int64},values["volume"])
        )
    
        DataFrames.rename!(new,[:Data, Symbol(symb)])
        #DataFrames.deleterows!(new, isnothing.(new).Close)
        if i == 1
            global old = new
        else
            old = outerjoin(old, new, on = :Data)
        end
    end
    return old
end
yahoo(s::Vector{Symbol}, date1, date2, interval) = yahoo(String.(s), date1, date2, interval)
yahoo(s::Vector{Symbol}, date1, date2) = yahoo(String.(s), date1, date2, "1d")
yahoo(s::Vector{Symbol}, date1) = yahoo(String.(s), date1, Date(Dates.now()), "1d")
yahoo(s::Symbol, date1, date2, interval) = yahoo(String(s), date1, date2, interval)
yahoo(s::Symbol, date1, date2) = yahoo(String(s), date1, date2, "1d")
yahoo(s::Symbol, date1) = yahoo(String(s), date1, Date(Dates.now()), "1d")


"""
    retornos(precos::DataFrame, tipo_retorno)
    
Calcular retornos (simples, ou logaritmicos) a partir de tabela de precos, com datas na primeira coluna

## Argumentos
* `precos`: tabela de precos, com datas na primeira coluna
* `tipo_retorno`: :simples, :log
"""
function retornos(df::DataFrame, tipo_retorno::Symbol) 
    c = Matrix(df)
    
    cot = convert(Matrix{Union{Missing, Float64}},c[:,2:end])
    dt  = convert(Vector{Date}, c[:,1])
    nret = size(c,1) - 1
    ncol = size(c,2) - 1       
    
    ret   = Array{Union{Missing, Float64}}(undef, nret, ncol)  
    datas = Array{Date}(undef,nret)
    if tipo_retorno == :log
        for i = 1:nret
            ret[i,:] = log.(cot[i+1,:]./cot[i,:])
            datas[i] = dt[i+1] 
        end
    elseif tipo_retorno == :simples
        for i = 1:nret
            ret[i,:] = cot[i+1,:]./cot[i,:] .- 1.0
            datas[i] = dt[i+1]    
        end
    end
    d = DataFrame(Data = datas)
    nomes = names(df)
    popfirst!(nomes)   
    r = DataFrame(ret, Symbol.(nomes))
    
    return  hcat(d, r)
end


function ETL(retornos, alfa, V)
    nobs = length(retornos)
    corte = quantile(retornos, alfa)
    excedimentos = 0
    n_excedimentos = 0
    for i in 1:nobs
        if retornos[i] <= corte
            excedimentos += retornos[i]
            n_excedimentos += 1
        end
    end
    return - (excedimentos/n_excedimentos - mean(retornos)) * V
end


function EWMA(retornos, lambda)
    m = length(retornos)
    media = mean(retornos)                 
    sigma2 = 0
    for k in 1:m
        sigma2 += (1-lambda)*lambda^(k-1) * (retornos[end-k+1] - media)^2
    end
    return sqrt(sigma2)                     
end


function semaforo(retornos, alfa)
    nobs = length(retornos)
    excep = zeros(nobs)
    sigma = std(retornos)
    RaR = - sigma * quantile(Normal(0,1), alfa)
    for i in 1:nobs
        if retornos[i] <= -RaR
            excep[i] = 1
        end
    end
    return sum(excep)
end


function cvm(Sigma)
    n = size(Sigma,1)
    modelo = Model(OSQP.Optimizer)
    @variable(modelo, w[1:n])                      
    @objective(modelo, Min, w' * Sigma * w) 
    @constraint(modelo, sum(w) == 1)        
    @constraint(modelo, w .>= 0)            
    optimize!(modelo)                       
    return value.(w)
end


function fe(mu, Sigma, mu_k)
    n = size(Sigma,1)
    modelo = Model(OSQP.Optimizer)
    @variable(modelo, w[1:n])     
    @objective(modelo, Min, w' * Sigma * w)  
    @constraint(modelo, sum(w) == 1)  
    @constraint(modelo, w .>= 0)      
    @constraint(modelo, mu' * w >= mu_k)
    optimize!(modelo)                   
    return value.(w)                                 
  end



function gfe(mu, Sigma, T, ncarteiras = 10)
    n = size(Sigma,1) 
    modelo = Model(OSQP.Optimizer)
    set_silent(modelo)                                 
    @variable(modelo, w[1:n])                           
    @objective(modelo, Min, w' * Sigma * w)             
    @constraint(modelo, sum(w) == 1)                    
    @constraint(modelo, w .>= 0)                        
    optimize!(modelo)                                   
    w_CVM = value.(w)  
    mu_CVM = mu' * w_CVM                                 
    var_k = zeros(ncarteiras)
    mu_k = zeros(ncarteiras)
    for i in 0:ncarteiras-1
        modelo = Model(OSQP.Optimizer)
        set_silent(modelo)                                  
        @variable(modelo, w[1:n])                           
        @objective(modelo, Min, w' * Sigma * w)             
        @constraint(modelo, sum(w) == 1)                    
        @constraint(modelo, w .>= 0)                        
        mu_max = maximum(mu)
        incremento = (mu_max-mu_CVM)/(ncarteiras-1)
        mu_k[i+1] = mu_CVM + incremento * i
        @constraint(modelo, mu' * w >= mu_k[i+1])                
        optimize!(modelo)                                  
        w_k = value.(w)   
        var_k[i+1] = T * w_k' * Sigma * w_k 
    end
    avar = T * diag(Sigma)
    fig = plot(var_k,mu_k, xlabel = "Risco da Carteira", ylabel = "Valor Esperado do Retorno da Carteira", label = "Fronteira Eficiente", xlim = (0, maximum(avar) * 1.1), ylim = (0, maximum(mu)*1.1), framestyle = :box, legend = :bottomright)
    fig = scatter!(avar, mu, label = "Ativos")
    var_CVM = T * w_CVM' * Sigma * w_CVM
    fig = scatter!([var_CVM], [mu_CVM], label = "Carteira variância minima")
    return fig
end

function alocar(mu, Sigma, lista, ncarteiras = 10)
    modelo = Model(OSQP.Optimizer)
    set_silent(modelo)       
    n = length(mu)                           
    @variable(modelo, w[1:n])                           
    @objective(modelo, Min, w' * Sigma * w)             
    @constraint(modelo, sum(w) == 1)                    
    @constraint(modelo, w .>= 0)                        
    optimize!(modelo)                                   
    w_CVM = value.(w)                                   
    mu_k = zeros(ncarteiras)
    w_k = zeros(n,ncarteiras)
    mu_CVM = mu' * w_CVM
    for i in 0:ncarteiras-1
        modelo = Model(OSQP.Optimizer)
        set_silent(modelo)                                  
        @variable(modelo, w[1:n])                           
        @objective(modelo, Min, w' * Sigma * w)             
        @constraint(modelo, sum(w) == 1)                    
        @constraint(modelo, w .>= 0)                        
        mu_max = maximum(mu)
        incremento = (mu_max-mu_CVM)/(ncarteiras-1)
        mu_k[i+1] = mu_CVM + incremento * i
        @constraint(modelo, mu' * w >= mu_k[i+1])                
        optimize!(modelo)                                  
        w_k[:,i+1] = value.(w)   
    end
    ticklabel = "P" .* string.(collect(1:ncarteiras))
    ticklabel[1] = "CVM"
    ticklabel[ncarteiras] = "CRM"
    fig = groupedbar(w_k', bar_position = :stack, bar_width=1.0, xlabel = "Carteiras", xticks=(1:ncarteiras, ticklabel), ylabel = "Pesos", label = permutedims(String.(lista)), legend = :bottomright, framestyle = :box)
    return fig
end


#=
function installedx()
    @warn "Pkg.installed() is deprecated"
    deps = dependencies()
    installs = Dict{String, VersionNumber}()
    for (uuid, dep) in deps
        dep.is_direct_dep || continue
        dep.version === nothing && continue
        installs[dep.name] = dep.version
    end
    return installs
end
=#
