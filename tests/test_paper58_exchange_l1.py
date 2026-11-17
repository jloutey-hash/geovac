"""Paper 58 -- the two-centre exchange integral with a 2p orbital on A (record EP, scored GO).

Guard for record docs/predictions/2026-11-17_exchange_l1.md (EP; owner Paper 58, row 3c of
docs/chemistry_period_inventory.md).  Written as a SEPARATE pass from the derivation and scoring (CLAUDE.md Sec. 9).

Two families, hydrogenic orbitals, A-orbital rate a, B-orbital (1s) rate b:

  S = (2p-sigma_A 1s_B | 2p-sigma_A 1s_B)     sigma = 0, prefactor a^5 b^3 R^7 / 32
  P = (2p_{+1,A} 1s_B | 2p_{-1,A} 1s_B)       sigma = 1, prefactor a^5 b^3 R^7 / 64

at (a, b) = (1, 1) and (1, 2).  Claims (the record's Results table):

  * operators over Q of order 4 / 4 / 9 / 9 and degree 16 / 9 / 34 / 25 (S11, P11, S12, P12) annihilate f;
  * R = 0 Fuchsian, exponents {-3, 0, 4, 11} at (1, 1), {-3 (x2), -1, 0, 1, 2, 4, 6, 11} at (1, 2);
  * rank one at infinity, rates {0, 2, -2 (x2)} at (1, 1), {-3, 4, 2, 0 (x2), -2 (x2), -4 (x2)} at (1, 2);
  * f's ln R part starts at R^11 ln R with coefficient 704, 16, 22528, 512 (/1157625); ratio (1,2)/(1,1) = 32.

PROVENANCE OF THE FROZEN OPERATORS.  `operator_Q_crt` of debug/data/ep_heldout_<case>_<a>_<b>_EPS.json (CRT +
rational reconstruction; sha256 2psig_1_1 22eb3c2a...891e, 2ppi_1_1 1582be53...1d49, 2psig_1_2 b85cd065...a6,
2ppi_1_2 52784199...7868), copied verbatim.  Nothing here reads or imports debug/.

INDEPENDENCE.  The exact annihilation tests use a Q-only port of the frozen record generator debug/ep_gen.py
(`generate`, record sha256 498496b3...7ecfa) together with the pieces it imports from debug/exch_gen.py; that is the
same series family the operators were fitted to (mod p), so it is NOT the independent check.  The independent checks
are (i) at (1, 1), where the tau sum terminates, the ODE residual L f = 0 evaluated on the production closed form
geovac.qfd_core.exchange_closed_form (symbolic in R, no series, no generator), and the generator's Frobenius value
against the same closed form; (ii) at (1, 2), the Frobenius solution of the frozen operator against the production
route geovac.qfd_core.exchange_hp (slow tier).

What each test rejects is named in its docstring; fire-test plants are listed at the end of the file.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import comb, factorial
from typing import Dict, List, Tuple

import numpy as np
import pytest
import sympy as sp

R = sp.Symbol("R")

# ---------------------------------------------------------------- frozen operators
# Coefficient lists [p_0, p_1, ..., p_r] of  L = sum_i p_i(R) D^i.

OP_S_1_1 = [
    "8*R*(2560*R**14 + 19200*R**13 + 116608*R**12 + 626672*R**11 + 2990512*R**10 + 10533652*R**9 + 26562666*R**8 + 44924520*R**7 + 996864*R**6 - 175063950*R**5 - 461929383*R**4 - 639167760*R**3 - 447650280*R**2 - 192099600*R - 96049800)",
    "-8*(512*R**16 + 2816*R**15 + 6528*R**14 + 12336*R**13 - 29936*R**12 - 927164*R**11 - 6609286*R**10 - 27789980*R**9 - 86941310*R**8 - 201967696*R**7 - 354039228*R**6 - 456853950*R**5 - 411941775*R**4 - 292347720*R**3 - 188862660*R**2 - 64033200*R - 32016600)",
    "-2*R*(1024*R**15 + 17920*R**14 + 111616*R**13 + 585312*R**12 + 2685968*R**11 + 9746232*R**10 + 26366200*R**9 + 48545204*R**8 + 67018156*R**7 + 76373106*R**6 + 92967210*R**5 + 132449205*R**4 + 146023920*R**3 + 115700760*R**2 + 64033200*R + 32016600)",
    "2*R**2*(512*R**14 + 256*R**13 - 11136*R**12 - 104400*R**11 - 683232*R**10 - 3414132*R**9 - 13379642*R**8 - 35776936*R**7 - 67904376*R**6 - 90126582*R**5 - 80089275*R**4 - 55954080*R**3 - 36514800*R**2 - 12806640*R - 6403320)",
    "R**3*(512*R**13 + 5120*R**12 + 34944*R**11 + 192624*R**10 + 857760*R**9 + 3019068*R**8 + 7753834*R**7 + 14622985*R**6 + 19833534*R**5 + 18957855*R**4 + 14700000*R**3 + 10195920*R**2 + 4268880*R + 2134440)",
]

OP_P_1_1 = [
    "24*R*(96*R**7 + 704*R**6 + 3352*R**5 + 12148*R**4 + 26180*R**3 + 36190*R**2 + 34650*R + 17325)",
    "-8*(96*R**9 + 608*R**8 + 1688*R**7 + 2204*R**6 - 7560*R**5 - 48300*R**4 - 131460*R**3 - 227430*R**2 - 242550*R - 121275)",
    "-2*R*(192*R**8 + 2656*R**7 + 14192*R**6 + 44240*R**5 + 101640*R**4 + 183960*R**3 + 253680*R**2 + 242550*R + 121275)",
    "2*R**2*(96*R**7 + 320*R**6 - 1288*R**5 - 9100*R**4 - 25620*R**3 - 45150*R**2 - 48510*R - 24255)",
    "R**3*(96*R**6 + 896*R**5 + 3640*R**4 + 9100*R**3 + 15330*R**2 + 16170*R + 8085)",
]

OP_S_1_2 = [
    "96*(445510172085608448000*R**32 + 12843396958531803955200*R**31 + 155211915887545203302400*R**30 + 1113783597300974769561600*R**29 + 5540221774795802630031360*R**28 + 20775255680503180898025984*R**27 + 61131331250646694307898624*R**26 + 140998779975431536511106816*R**25 + 241376164933680775886908032*R**24 + 245955015511884705501456768*R**23 - 79433595218807165475746496*R**22 - 888975555506667903623172288*R**21 - 1808226002182999702063117824*R**20 - 1573805929599509083847597248*R**19 + 1263127415103927652253663808*R**18 + 6653857662987335897513537064*R**17 + 11990264904013217317844500656*R**16 + 13374575478088562991857588364*R**15 + 8912781385135025269864396824*R**14 + 603637578231467858155811830*R**13 - 6871088286729935584765129980*R**12 - 9659421429094444912952144889*R**11 - 7311434283968457581948098908*R**10 - 2614525150668493057222224186*R**9 + 911628005494031238436503696*R**8 + 1968890396033303685177875472*R**7 + 1441789729378972965042134760*R**6 + 711510833496561582910565340*R**5 + 309124878174351620007777000*R**4 + 120300197542370409301531200*R**3 + 33098714817865332408907200*R**2 + 5016524096142995635584000*R + 282723797213228370384000)",
    "-96*(207904746973283942400*R**33 + 5736179370109823631360*R**32 + 64371637294655529000960*R**31 + 413114881417575746273280*R**30 + 1752625220205263918376960*R**29 + 5225012992261184480745984*R**28 + 10546310043060585793691904*R**27 + 9113046428170321530326784*R**26 - 29110311703937544277756032*R**25 - 161186144114047756365076608*R**24 - 423615285179952546180627648*R**23 - 680733398501847196268620992*R**22 - 464237531551582140174395328*R**21 + 938765727171996856071449152*R**20 + 3765479964588725628037261824*R**19 + 6978887620534144524704153576*R**18 + 8465886489648420986030371872*R**17 + 6500063871634644872832832940*R**16 + 1599718134288286214573016600*R**15 - 3492552102034396730856063946*R**14 - 5900205196623193480224340596*R**13 - 4852519041458474514498127817*R**12 - 2041342940457306763537472208*R**11 + 262724376080515656770587510*R**10 + 1095111288993458567376318600*R**9 + 950753709642128686622779536*R**8 + 571520918441395578686529384*R**7 + 297461116874275473421715388*R**6 + 140675015723617244191630680*R**5 + 48896026432459954454772960*R**4 + 7772938894979085812149440*R**3 - 1404537330213841510874880*R**2 - 895293171443325332937600*R - 128383480164916040486400)",
    "24*(118802712556162252800*R**34 + 2781391162595935518720*R**33 + 21161261408064091914240*R**32 + 29083503952701513768960*R**31 - 606068217997206627962880*R**30 - 5534777786645488507849728*R**29 - 27581148425717526055735296*R**28 - 97796633443609797283692672*R**27 - 263620464340884994642731840*R**26 - 539320668108225728478721344*R**25 - 768490821394119645959568096*R**24 - 448453392414159964335486144*R**23 + 1190271412145452771834883280*R**22 + 4523154235704899145701371240*R**21 + 8353687113125459250913579248*R**20 + 9695196197833170670707972692*R**19 + 5780380987187775074771853984*R**18 - 3081126362206826672544287560*R**17 - 12288034662921666970384615632*R**16 - 15795591834438479953791002968*R**15 - 11165261489920150133541529104*R**14 - 1831888002852295130975537780*R**13 + 5530221485116012100638130796*R**12 + 7080153744647513485517003143*R**11 + 4393196787320223243945070356*R**10 + 1232803988201924195653336686*R**9 - 326610350232771755765293920*R**8 - 595138298505552191832680184*R**7 - 425035600161427918793667240*R**6 - 233100729906521557915426260*R**5 - 112143281758133120593811640*R**4 - 44689176230265249215037120*R**3 - 12135965236236217619496000*R**2 - 1724419047314870005017600*R - 68243675189399951472000)",
    "8*(386108815807527321600*R**34 + 12123441691873835581440*R**33 + 157926121178162779054080*R**32 + 1177593005759575408680960*R**31 + 5819760361169476585113600*R**30 + 20613718718814148806979584*R**29 + 53233784080621430294717568*R**28 + 92013779541198260951106624*R**27 + 51879184754935825047320640*R**26 - 295036470935520354056925504*R**25 - 1285911194485971868125404928*R**24 - 3008671236655319718348235488*R**23 - 4607998442023495591408663288*R**22 - 3903941289368022020835829272*R**21 + 1483049186535613258845232692*R**20 + 11442976108700184729476259372*R**19 + 21805068807634630164153783840*R**18 + 26175372439577175029275041762*R**17 + 20547868415224699637856320332*R**16 + 7145392549981386240389304495*R**15 - 6560815762149825806523228048*R**14 - 13127282431781606112736781424*R**13 - 11147040453293245855921344189*R**12 - 5201857612329803647663480560*R**11 - 385893386776437739497999102*R**10 + 1454009607661592849258826360*R**9 + 1482677730734937421905343656*R**8 + 1003437095026489623340131096*R**7 + 530988113278988890820117532*R**6 + 227672868844343392710837480*R**5 + 76941710882111590506616320*R**4 + 16715279424800597450420160*R**3 + 931281904868833402663680*R**2 - 436840837701906221654400*R - 71285762585567341785600)",
    "12*R*(14850339069520281600*R**33 + 1186718052752387850240*R**32 + 26201699346991793080320*R**31 + 280655883236500408189440*R**30 + 1868786663183720612663040*R**29 + 8845509362700143232918144*R**28 + 31941965535784121023831872*R**27 + 90892280766889678851213120*R**26 + 204723878911420253553866016*R**25 + 353700344708066480294444736*R**24 + 415015214804170190091570672*R**23 + 130485252376062207320817736*R**22 - 772907640601212167629881888*R**21 - 2277133276058984115928588164*R**20 - 3777698879290188943081225272*R**19 - 4265532253226498332278481564*R**18 - 3080416667227669953851279940*R**17 - 541790314594875227701239386*R**16 + 2097872134462798678734587946*R**15 + 3442407617488639609589608475*R**14 + 2957257089863323470847605096*R**13 + 1382761452617203017798533949*R**12 - 12285955591731745087501764*R**11 - 559021947752708969253118986*R**10 - 466947485311018653453260808*R**9 - 268799617561585617337572636*R**8 - 162902105882652514394343744*R**7 - 93429979998772528468777428*R**6 - 39248503905515442572163960*R**5 - 10008373925583540104904000*R**4 - 295519874626357248162240*R**3 + 932386470515759682593280*R**2 + 359470900324004099241600*R + 47523841723711561190400)",
    "-12*R**2*(70539110580221337600*R**32 + 1938866940199343831040*R**31 + 21664062831864892331520*R**30 + 140206917460344020709120*R**29 + 611249664417348261018240*R**28 + 1902038012108666343261504*R**27 + 4104945678174031545758784*R**26 + 4533272392185131789907648*R**25 - 7333289356783129114239360*R**24 - 54132377544386960197315968*R**23 - 163308242387955951873553128*R**22 - 331314453018673111963596312*R**21 - 472143565627897682958530892*R**20 - 411231158759608187135609988*R**19 - 4010296033454665589107076*R**18 + 673817235740727020307710370*R**17 + 1302080912347441330177306298*R**16 + 1518422998074412538623086711*R**15 + 1179147810770404558337306491*R**14 + 492872360370436550996998344*R**13 - 129752650249282903652649997*R**12 - 398058220073543445270662976*R**11 - 331375556426164295720897850*R**10 - 153695760072614693101010304*R**9 - 41612483358679514816517372*R**8 - 8420815502835199074763488*R**7 - 426381005091038926790940*R**6 + 3493702393684301017857816*R**5 + 3257955836970156999512064*R**4 + 1409383945591077144035520*R**3 + 365808546673447559489280*R**2 + 67949405131863344158080*R + 7920640287285260198400)",
    "-3*R**3*(96527203951881830400*R**31 + 3011678735003661864960*R**30 + 39780359201128877061120*R**29 + 312810815785589912240640*R**28 + 1701395108290635021515520*R**27 + 6934828342066823285006976*R**26 + 22099478060715601263155712*R**25 + 55914641909182094296124160*R**24 + 111076321673743384488920640*R**23 + 164088133777419901391083488*R**22 + 145122874202647502171661600*R**21 - 42022638193144123452257272*R**20 - 431294826413667742321855968*R**19 - 886648800202366757209719292*R**18 - 1110477091448692921381285536*R**17 - 862819132441426427960166002*R**16 - 200346866689164684576114060*R**15 + 539402463941725302292473084*R**14 + 978942491982752427690445116*R**13 + 972136122996452585700648347*R**12 + 656317030656294592106977764*R**11 + 281335900117195371287202542*R**10 + 25612919370308605705612272*R**9 - 80108429433945747003251208*R**8 - 90240568178808123584989416*R**7 - 63461272286759644942892772*R**6 - 32791628214435699741527064*R**5 - 11401054323615318143752896*R**4 - 1763330534946674365273920*R**3 + 386040360461704476549120*R**2 + 230244804967685406180480*R + 31682561149141040793600)",
    "3*R**4*(3712584767380070400*R**30 + 58145941883927439360*R**29 - 105286468114229921280*R**28 - 7484323443843150408960*R**27 - 73046589805067261097600*R**26 - 414448620692873564119104*R**25 - 1668644900418211204128768*R**24 - 5152275807358760645992704*R**23 - 12612546164784885601647168*R**22 - 24593036269848050092604112*R**21 - 37322327709260885263814664*R**20 - 40444510792680278372408832*R**19 - 20201041922601686147790060*R**18 + 28720428533219506003890588*R**17 + 89139096096724716905267262*R**16 + 125429774362096077161282664*R**15 + 111302577795033170394852172*R**14 + 54133939631789157333371670*R**13 - 8146214824619700247771721*R**12 - 39951039206742028242107832*R**11 - 36356822333848681052796254*R**10 - 17637370344274490329271004*R**9 - 4992440960805840826264032*R**8 - 1743246049147881692717736*R**7 - 1212193417930236500004468*R**6 - 516925617579512861272632*R**5 - 37804029108248844307008*R**4 + 75764388654973638244800*R**3 + 47198878909803318551040*R**2 + 13920347495165188429440*R + 1760142286063391155200)",
    "3*R**5*(5568877151070105600*R**29 + 164255046749027619840*R**28 + 2019871132086683646720*R**27 + 14687408147994174366720*R**26 + 73763448587638581985920*R**25 + 277379564125532403749184*R**24 + 812146582956586205637072*R**23 + 1866197611733317579659600*R**22 + 3273093862419453580321200*R**21 + 3962030545626324707779032*R**20 + 1847501216871682285759140*R**19 - 4699961919848295287375432*R**18 - 14047835317131023725795344*R**17 - 20496449755270969751278234*R**16 - 17903562468617097084797649*R**15 - 5307922470457080507240089*R**14 + 10780351540148441867458347*R**13 + 20936685306607300348827932*R**12 + 20503078559532306522343251*R**11 + 12650483893377908708990985*R**10 + 4237126249245821146726236*R**9 - 523210208104866883762620*R**8 - 1893904713244714153240836*R**7 - 1688736920215939939936704*R**6 - 1008445403404102002180588*R**5 - 385335729137355834093900*R**4 - 65326087223416125278640*R**3 + 11681863439687395540080*R**2 + 7585929002767247779680*R + 1045084482350138498400)",
    "R**6*(1856292383690035200*R**28 + 56607974633365908480*R**27 + 723608577036330912000*R**26 + 5458609387026315498240*R**25 + 28226889120563737660800*R**24 + 108141459775490655505728*R**23 + 318778398663524582103792*R**22 + 728325336798947387254464*R**21 + 1252881362317361723941392*R**20 + 1459885888799592871919832*R**19 + 615833738957227428586380*R**18 - 1635079944278012587856964*R**17 - 4319260673648560666852456*R**16 - 5392396360540803028141926*R**15 - 3571233551520897682425027*R**14 + 214711149581471876711652*R**13 + 3450309746995166038011681*R**12 + 4295320854539526531155781*R**11 + 3016296538031189701975945*R**10 + 1200531426428113201021050*R**9 + 78451609272925714440912*R**8 - 270288702924087056681988*R**7 - 270916633849455308837628*R**6 - 171588280067211200546772*R**5 - 69233365474594355431884*R**4 - 13056300588128018978160*R**3 + 1435838637755547534960*R**2 + 1181814830801986502880*R + 165013339318442920800)",
]

OP_P_1_2 = [
    "384*(1325049062400*R**23 + 63025012189440*R**22 + 984457696320000*R**21 + 8741094206872320*R**20 + 55472607026555040*R**19 + 268455546394575624*R**18 + 1020769860414328416*R**17 + 3114127071498931920*R**16 + 7627266270985844520*R**15 + 14965002045766458468*R**14 + 23398402329817794082*R**13 + 28502192458888555752*R**12 + 25240415228194056128*R**11 + 12278605061678931972*R**10 - 5244337287523359122*R**9 - 18343667024712689202*R**8 - 21433920463255795633*R**7 - 16593076945410460110*R**6 - 9690315481686017690*R**5 - 4487214477176953200*R**4 - 1583837706173510700*R**3 - 358850828220687000*R**2 - 34510674922911000*R + 1294460788104000)",
    "-64*(3710137374720*R**24 + 171876530714112*R**23 + 2517380405987328*R**22 + 20138010468629760*R**21 + 111659134845003840*R**20 + 453605338879313328*R**19 + 1341184673987855712*R**18 + 2708445939201990576*R**17 + 2440441822017536784*R**16 - 5667386493201229296*R**15 - 30948945602402824224*R**14 - 78951328836069555414*R**13 - 140479195145958571816*R**12 - 188275079952284681088*R**11 - 191322319200274880112*R**10 - 139030576238268321642*R**9 - 52685812508954693568*R**8 + 30395499241032876735*R**7 + 81418045587262696590*R**6 + 90748916691996621150*R**5 + 68448860141839115400*R**4 + 36946755832575827700*R**3 + 14040669325294438200*R**2 + 3409869867614397000*R + 399464850000744000)",
    "16*(2120078499840*R**25 + 89356260962304*R**24 + 1001455020195840*R**23 + 4321988544570624*R**22 - 1171337653873152*R**21 - 143407940265721440*R**20 - 1100703460660531056*R**19 - 5129229547979372412*R**18 - 17276755481510156568*R**17 - 43931483338481848800*R**16 - 85904399646948191244*R**15 - 128704223340526989564*R**14 - 139223593229165449918*R**13 - 81989342383679137800*R**12 + 49772139155681573880*R**11 + 215118399976260267060*R**10 + 344318168969711312274*R**9 + 386609551270560048330*R**8 + 339470121127618625067*R**7 + 240083109772698523410*R**6 + 137219961475693747470*R**5 + 63415357799948625600*R**4 + 23316213570608960100*R**3 + 6088078548365913000*R**2 + 873135226621845000*R + 27183676550184000)",
    "16*(2296751708160*R**25 + 111260373590016*R**24 + 1757576089979136*R**23 + 14997666727423872*R**22 + 86871514717441152*R**21 + 373701911812987968*R**20 + 1216607854519378620*R**19 + 2997611200907378964*R**18 + 5387344627704228918*R**17 + 6148341702400346972*R**16 + 675375874335394524*R**15 - 16294854243559947774*R**14 - 46128497201088814782*R**13 - 81994668230405208198*R**12 - 109156721993602109775*R**11 - 112094500214975817968*R**10 - 85054568093199964692*R**9 - 38034515831571689527*R**8 + 7107007004022823326*R**7 + 31841497210188490314*R**6 + 33883359423091113960*R**5 + 23702528225324521260*R**4 + 12223378291662232200*R**3 + 4540786597022072760*R**2 + 1083542033792433600*R + 123413967815846400)",
    "8*R*(265009812480*R**24 + 19252331900928*R**23 + 486833775301440*R**22 + 5839649997715440*R**21 + 42457631315085360*R**20 + 221249484058794084*R**19 + 891295872647313600*R**18 + 2885510462681201148*R**17 + 7753534564025878416*R**16 + 17608655473769896524*R**15 + 33932544081076808682*R**14 + 55037376383254891308*R**13 + 73664516208719099704*R**12 + 77873436099672172182*R**11 + 57631123865281749680*R**10 + 13910558329862636814*R**9 - 38078274963823955233*R**8 - 76990554358288539462*R**7 - 89612920669525185258*R**6 - 77299430392497636720*R**5 - 51491998831010013420*R**4 - 25990930896447389400*R**3 - 9422088620664998520*R**2 - 2194267744135051200*R - 246827935631692800)",
    "-8*R**2*(1258796609280*R**23 + 58758511952448*R**22 + 863849677621392*R**21 + 6863603058522432*R**20 + 38371153590336348*R**19 + 165355915454936172*R**18 + 551387461206405258*R**17 + 1417467517542546156*R**16 + 2701437028931649870*R**15 + 3295696493508258810*R**14 + 454221601780069278*R**13 - 9415311749405272314*R**12 - 27927720392032336347*R**11 - 51500827585925377932*R**10 - 71343984388940808309*R**9 - 78411941804988022293*R**8 - 69938080097177842596*R**7 - 51266354296171345718*R**6 - 31040731055539967364*R**5 - 15342023121722726316*R**4 - 5994901724518201680*R**3 - 1777976077837239720*R**2 - 368912923705085760*R - 41137989271948800)",
    "-2*R**3*(1722563781120*R**22 + 82705461132672*R**21 + 1275622832423520*R**20 + 10660356355607712*R**19 + 61983398363526936*R**18 + 275340131644652928*R**17 + 947797471471175136*R**16 + 2567273405244464070*R**15 + 5490711052458524442*R**14 + 9057095990343844308*R**13 + 10653668049852968250*R**12 + 6034604045259411198*R**11 - 8148398867898934366*R**10 - 30741379861116673854*R**9 - 54631699173401678409*R**8 - 69807696354144867582*R**7 - 69680015980916210842*R**6 - 55579820169031932384*R**5 - 35315488505201084076*R**4 - 17403553188832465080*R**3 - 6261117763386392280*R**2 - 1459099856209786560*R - 164551957087795200)",
    "2*R**4*(66252453120*R**21 + 2897282872512*R**20 + 37734169343472*R**19 + 249767222914224*R**18 + 1162056258314316*R**17 + 4780536848136504*R**16 + 16769820589156566*R**15 + 38884302832108914*R**14 + 32138857545892992*R**13 - 129026676348437226*R**12 - 636769097722976976*R**11 - 1599754179106432374*R**10 - 2853754280898834042*R**9 - 3953730437832202593*R**8 - 4433807564947171458*R**7 - 4127631633384580754*R**6 - 3209250587031004752*R**5 - 2031806782424398188*R**4 - 996082207735704840*R**3 - 354167041114343160*R**2 - 81510808389151680*R - 9141775393766400)",
    "R**5*(198757359360*R**20 + 9387499375296*R**19 + 139558732926480*R**18 + 1099010017622520*R**17 + 5984772235386156*R**16 + 25232176114935246*R**15 + 82065904110243918*R**14 + 201211343115607224*R**13 + 360569571982322184*R**12 + 425211182948453922*R**11 + 153326166062273784*R**10 - 640108690121588301*R**9 - 1838524737214036083*R**8 - 2956956247737760158*R**7 - 3430142496593341726*R**6 - 3047150500436451708*R**5 - 2095494516776379060*R**4 - 1091894525388239400*R**3 - 406621645645777320*R**2 - 96252928438201920*R - 10855858280097600)",
    "R**6*(22084151040*R**19 + 1043055486144*R**18 + 15390630826704*R**17 + 118692083996568*R**16 + 625410664821828*R**15 + 2525614828405326*R**14 + 7715314440913032*R**13 + 17213272941125448*R**12 + 26675184599382672*R**11 + 23534411794821394*R**10 - 6498171121235418*R**9 - 63902997656581569*R**8 - 126176862554626536*R**7 - 160314877453582414*R**6 - 149561009966307816*R**5 - 105921817878682476*R**4 - 56296360955126880*R**3 - 21239194456245480*R**2 - 5061704432178240*R - 571360962110400)",
]

# key: (family, a, b, operator, order, degree)
CASES_T = {
    "S_1_1": ("2psig", Fraction(1), Fraction(1), OP_S_1_1, 4, 16),
    "P_1_1": ("2ppi", Fraction(1), Fraction(1), OP_P_1_1, 4, 9),
    "S_1_2": ("2psig", Fraction(1), Fraction(2), OP_S_1_2, 9, 34),
    "P_1_2": ("2ppi", Fraction(1), Fraction(2), OP_P_1_2, 9, 25),
}
KEYS = list(CASES_T)
KEYS_EQ = ["S_1_1", "P_1_1"]
KEYS_UNEQ = ["S_1_2", "P_1_2"]

EXPONENTS_EQ = {-3: 1, 0: 1, 4: 1, 11: 1}
EXPONENTS_UNEQ = {-3: 2, -1: 1, 0: 1, 1: 1, 2: 1, 4: 1, 6: 1, 11: 1}
LN_COEF = {"S_1_1": Fraction(704, 1157625), "P_1_1": Fraction(16, 1157625),
           "S_1_2": Fraction(22528, 1157625), "P_1_2": Fraction(512, 1157625)}


def exponents0(key: str) -> Dict[int, int]:
    return EXPONENTS_EQ if CASES_T[key][1] == CASES_T[key][2] else EXPONENTS_UNEQ


def rates_inf(a: Fraction, b: Fraction) -> Dict[sp.Rational, int]:
    """Equal exponents: {0, 2a, -2a (x2)}; unequal: {-(a+b), 2a, 2b, 0 (x2), -2a (x2), -2b (x2)}."""
    q = lambda x: sp.Rational(x.numerator, x.denominator)  # noqa: E731
    if a == b:
        return {sp.Integer(0): 1, q(2 * a): 1, q(-2 * a): 2}
    return {q(-(a + b)): 1, q(2 * a): 1, q(2 * b): 1, sp.Integer(0): 2, q(-2 * a): 2, q(-2 * b): 2}


@lru_cache(maxsize=None)
def op_coeffs(key: str) -> Tuple[Tuple[Fraction, ...], ...]:
    """p_i as ascending integer coefficient tuples."""
    out = []
    for s in CASES_T[key][3]:
        P = sp.Poly(sp.sympify(s, locals={"R": R}), R)
        cs = [Fraction(0)] * (P.degree() + 1)
        for (j,), c in P.terms():
            cs[j] = Fraction(int(c.p), int(c.q))
        out.append(tuple(cs))
    return tuple(out)


# ---------------------------------------------------------------- series generator
# Q-only port of debug/ep_gen.py `generate` (record EP, frozen) and the debug/exch_gen.py pieces it imports
# (RingQ, Tables, _biv_eval_uminus1, _biv_eval_uhalf, _shift2, legendre_PV, leg_moment, _prime_logs).  Code kept
# ring-generic as in the source, with only the Q ring present.


class RingQ:
    name = "Q"
    p = None

    def vec(self, n: int) -> np.ndarray:
        a = np.empty(n, dtype=object)
        a[:] = [Fraction(0)] * n
        return a

    def mat(self, n: int, m: int) -> np.ndarray:
        a = np.empty((n, m), dtype=object)
        for i in range(n):
            a[i, :] = [Fraction(0)] * m
        return a

    def c(self, x) -> Fraction:
        if isinstance(x, np.integer):
            x = int(x)
        x = Fraction(x)
        return Fraction(int(x.numerator), int(x.denominator))

    def red(self, a):
        return a

    def conv(self, a, b):
        return np.convolve(a, b)

    def dot(self, a, B):
        return a.dot(B)


def legendre_PV(tmax: int):
    """P_n and V_n (Q_n = P_n Q_0 - V_n) coefficient lists, n = 0..tmax."""
    P = [[Fraction(1)], [Fraction(0), Fraction(1)]]
    V = [[Fraction(0)], [Fraction(1)]]
    for n in range(1, tmax):
        for L in (P, V):
            a, b = L[n], L[n - 1]
            new = [Fraction(0)] * (len(a) + 1)
            for i, c in enumerate(a):
                new[i + 1] += Fraction(2 * n + 1, n + 1) * c
            for i, c in enumerate(b):
                new[i] -= Fraction(n, n + 1) * c
            L.append(new)
    return P[: tmax + 1], V[: tmax + 1]


def leg_moment(n: int, tau: int) -> Fraction:
    """int_{-1}^1 x^n P_tau(x) dx."""
    if n < tau or (n - tau) % 2:
        return Fraction(0)
    return Fraction(2 ** (tau + 1) * factorial(n) * factorial((n + tau) // 2),
                    factorial((n - tau) // 2) * factorial(n + tau + 1))


def _stack_obj(rows):
    M = np.empty((len(rows), len(rows[0])), dtype=object)
    for i, v in enumerate(rows):
        M[i, :] = v
    return M


class Tables:
    """int_0^oo e^{-pv} v^n B(v) dv as Laurent arrays in p, n = 0..nmax, exponents emin..emax."""

    def __init__(self, ring, nmax: int, emax: int):
        self.ring, self.nmax, self.emax = ring, nmax, emax
        self.emin = -(nmax + 1)
        W = emax - self.emin + 1
        self.W = W
        r = ring
        fact = [factorial(n) for n in range(nmax + 1)]
        H = [Fraction(0)]
        for n in range(1, nmax + 1):
            H.append(H[-1] + Fraction(1, n))
        powT = r.mat(nmax + 1, W)
        lnv1 = r.mat(nmax + 1, W)
        for n in range(nmax + 1):
            i = -n - 1 - self.emin
            powT[n, i] = r.c(fact[n])
            lnv1[n, i] = r.c(fact[n] * H[n])
        neg = r.red(-powT)
        T = {"pow": {"1": powT}, "lnv": {"1": lnv1, "g": neg, "lp": neg.copy()}, "l2": {"l2": powT.copy()}}
        for a, lmult, key in ((2, 1, "a2"), (4, 2, "a4")):
            T[key] = self._tn(a, lmult)
        self.T = T

    def _tn(self, a: int, lmult: int):
        """T_n^a = int e^{-pv} v^n ln(v+a); T_0 = (1/p)[ln a + e^{ap} E1(ap)]; T_n = (-d/dp)^n T_0."""
        r = self.ring
        nmax, emin = self.nmax, self.emin
        top = self.emax + nmax + 2
        L = top + 2
        ea = [Fraction(a) ** k / factorial(k) for k in range(L + 1)]
        S = [Fraction(0)] + [Fraction(-a) ** k / (k * factorial(k)) for k in range(1, L + 1)]
        eS = [sum(ea[i] * S[k - i] for i in range(k + 1)) for k in range(L + 1)]
        Wx = top - emin + 1

        def arr(series_nonneg, shift):
            v = r.vec(Wx)
            for k, c in enumerate(series_nonneg):
                e = k + shift
                if emin <= e <= top and c != 0:
                    v[e - emin] = r.c(c)
            return v
        one = arr([-c for c in eS], -1)
        g = arr([-c for c in ea], -1)
        la = arr([Fraction(1) - ea[0]] + [-c for c in ea[1:]], -1)
        lp = arr([-c for c in ea], -1)
        em = np.array([Fraction(int(e)) for e in range(emin, top + 1)], dtype=object)

        def d(v):
            w = r.red(v * em)
            out = r.vec(Wx)
            out[:-1] = w[1:]
            return out

        rows = {k: [] for k in ("1", "g", "l2", "lp")}
        cur = {"1": one, "g": g, "la": la, "lp": lp}
        for _n in range(nmax + 1):
            for k in ("1", "g", "lp"):
                rows[k].append(cur[k][: self.W])
            rows["l2"].append(r.red(cur["la"][: self.W] * lmult))
            lp_over_p = r.vec(Wx)
            lp_over_p[:-1] = cur["lp"][1:]
            cur = {"1": r.red(-(d(cur["1"]) + lp_over_p)), "g": r.red(-d(cur["g"])),
                   "la": r.red(-d(cur["la"])), "lp": r.red(-d(cur["lp"]))}
        return {k: _stack_obj(v) for k, v in rows.items()}


def _biv_eval_uminus1(r, Tm: np.ndarray) -> np.ndarray:
    """sum_{i,k} T[i,k] u^i (u-1)^k -> coefficients in u."""
    Du, Dy = Tm.shape
    res = r.vec(Du + Dy)
    for k in range(Dy - 1, -1, -1):
        sh = r.vec(Du + Dy)
        sh[1:] = res[:-1]
        res = r.red(sh - res)
        res[:Du] = r.red(res[:Du] + Tm[:, k])
    return res


def _biv_eval_uhalf(r, Tm: np.ndarray) -> np.ndarray:
    """sum_{i,k} T[i,k] u^i (u/2)^k -> coefficients in u."""
    Du, Dy = Tm.shape
    res = r.vec(Du + Dy)
    ih = r.c(Fraction(1, 2))
    w = r.c(1)
    for k in range(Dy):
        res[k:k + Du] = r.red(res[k:k + Du] + r.red(Tm[:, k] * w))
        w = w * ih
    return res


def _shift2(r, c: np.ndarray) -> np.ndarray:
    """q(u) -> q(2+v) coefficients in v."""
    D = len(c)
    res = r.vec(D)
    two = r.c(2)
    for i in range(D - 1, -1, -1):
        sh = r.vec(D)
        sh[1:] = res[:-1]
        res = r.red(sh + r.red(res * two))
        res[0] = res[0] + c[i]
    return res


def _prime_logs(x: Fraction) -> Dict[int, int]:
    out: Dict[int, int] = {}
    for pr, e in sp.factorint(x.numerator).items():
        out[int(pr)] = out.get(int(pr), 0) + int(e)
    for pr, e in sp.factorint(x.denominator).items():
        out[int(pr)] = out.get(int(pr), 0) - int(e)
    return {k: v for k, v in out.items() if v}


GEN_CASES = {
    "2psig": {"sigma": 0, "H": 0, "c": {(2, 0): 1, (0, 2): -1, (3, 1): 1, (1, 3): -1}, "rexp": 7,
              "pref": lambda a, b: a ** 5 * b ** 3 / 32},
    "2ppi": {"sigma": 1, "H": 1, "c": {(2, 0): 1, (0, 2): -1}, "rexp": 7,
             "pref": lambda a, b: a ** 5 * b ** 3 / 64},
}


def pmul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def pder(a):
    return [i * a[i] for i in range(1, len(a))] or [Fraction(0)]


def ppow(a, n):
    out = [Fraction(1)]
    for _ in range(n):
        out = pmul(out, a)
    return out


def xpow(j):
    return [Fraction(0)] * j + [Fraction(1)]


SQ_M1 = [Fraction(-1), Fraction(0), Fraction(1)]       # x^2 - 1


def tau_polys(tau, sigma, H, js, Pc, Vc):
    """F_j, B_j, C_j for every j in js (G_j = B_j Q0 - C_j)."""
    P, V = [Fraction(x) for x in Pc[tau]], [Fraction(x) for x in Vc[tau]]
    if sigma == 0:
        DP, Cbase = P, pmul(ppow(SQ_M1, H), V)
    elif sigma == 1:
        DP = pder(P)
        Cbase = padd(pmul(ppow(SQ_M1, H), pder(V)), pmul(ppow(SQ_M1, H - 1), P))
    else:
        raise ValueError("sigma <= 1 only")
    base = pmul(ppow(SQ_M1, H), DP)
    F = {j: pmul(xpow(j), base) for j in js}
    return F, F, {j: pmul(xpow(j), Cbase) for j in js}


def eta_moment(n, tau, sigma, H):
    """int_{-1}^1 eta^n (1-eta^2)^H D^sigma P_tau(eta) d eta, exact."""
    if sigma == 0 and H == 0:
        return leg_moment(n, tau)
    if sigma == 1 and H == 1:
        out = Fraction(0)
        if n >= 1:
            out -= n * leg_moment(n - 1, tau)
        out += (n + 2) * leg_moment(n + 1, tau)
        return out
    raise ValueError("(sigma, H) not supported")


def _bivariate(r, Fc):
    D = len(Fc)
    A = r.mat(D, D)
    for m, cm in enumerate(Fc):
        if cm == 0:
            continue
        for i in range(m + 1):
            A[i, m - i] = A[i, m - i] + r.c(cm * comb(m, i) * (-1) ** (m - i))
    return A


def g_pieces_gen(r, Fd, Bd, Cd, j1, j2):
    """g(2+v) = int_{u/2}^{u-1} [F_j1(u-y) G_j2(y) + F_j2(u-y) G_j1(y)] dy, as polynomials in v times
    {1, ln v, ln(v+2), ln(v+4), ln 2}."""
    Du = max(len(Fd[j1]), len(Fd[j2]))
    LB = max(max(len(Bd[j]), len(Cd[j])) for j in (j1, j2))
    Dy2 = Du + LB - 1
    M = r.mat(Du, Dy2)
    N = r.mat(Du, Dy2)
    for ja, jb in ((j1, j2), (j2, j1)):
        K = _bivariate(r, Fd[ja])
        bv = np.array([r.c(x) for x in Bd[jb]], dtype=K.dtype)
        cv = np.array([r.c(x) for x in Cd[jb]], dtype=K.dtype)
        for i in range(K.shape[0]):
            mb = r.conv(K[i, :], bv)
            nb = r.conv(K[i, :], cv)
            M[i, :len(mb)] = M[i, :len(mb)] + mb
            N[i, :len(nb)] = N[i, :len(nb)] + nb
    Dy3 = Dy2 + 1
    ks = np.arange(Dy2)
    invk1 = np.array([r.c(Fraction(1, k + 1)) for k in range(Dy2)], dtype=M.dtype)
    s = M * invk1[None, :]
    Nint = r.mat(Du, Dy3)
    Nint[:, 1:] = N * invk1[None, :]
    out_lc, out_rc = {}, {}
    for c in (+1, -1):
        mc = -c
        L = r.mat(Du, Dy3)
        L[:, 1:] = s
        sign_k1 = np.array([r.c(mc ** (int(k) + 1)) for k in ks], dtype=M.dtype)
        L[:, 0] = L[:, 0] - (s * sign_k1[None, :]).sum(axis=1)
        w = s * sign_k1[None, :]
        suf = np.flip(np.cumsum(np.flip(w, axis=1), axis=1), axis=1)
        Rm = r.mat(Du, Dy3)
        for j in range(1, Dy3):
            Rm[:, j] = suf[:, j - 1] * r.c(Fraction(-1, j) * (mc ** j))
        out_lc[c], out_rc[c] = L, Rm
    half = r.c(Fraction(1, 2))
    src = (("L+", out_lc[1]), ("L-", out_lc[-1]), ("R+", out_rc[1]), ("R-", out_rc[-1]), ("N", Nint))
    up = {k: _biv_eval_uminus1(r, v) for k, v in src}
    lo = {k: _biv_eval_uhalf(r, v) for k, v in src}
    h = lambda a: a * half  # noqa: E731
    pieces = {
        "a2": h(up["L+"]),
        "a4": -h(lo["L+"]),
        "lnv": h(lo["L-"]) - h(up["L-"]),
        "l2": h(lo["L+"]) - h(lo["L-"]),
        "pow": h(up["R+"] - lo["R+"]) - h(up["R-"] - lo["R-"]) - (up["N"] - lo["N"]),
    }
    out = {}
    for k, v in pieces.items():
        nz = np.nonzero(v != 0)[0]
        top = int(nz[-1]) + 1 if len(nz) else 1
        out[k] = _shift2(r, v[:top])
    return out


def _vb(k, tau, sigma, H):
    return max(0, tau + sigma - k - 2 * H)


def _pair_terms(c):
    out: Dict[Tuple[int, int], List] = {}
    for (j1, k1), c1 in c.items():
        for (j2, k2), c2 in c.items():
            key = (min(j1, j2), max(j1, j2))
            out.setdefault(key, []).append((k1, k2, Fraction(c1 * c2)) if j1 <= j2 else (k2, k1, Fraction(c1 * c2)))
    return out


def tau_bound(case, tau):
    C = GEN_CASES[case]
    s, H = C["sigma"], C["H"]
    best = None
    for (j1, j2), lst in _pair_terms(C["c"]).items():
        vmin = -(j1 + j2 + 4 * H - 2 * s + 1)
        for k1, k2, _ in lst:
            v = C["rexp"] + _vb(k1, tau, s, H) + _vb(k2, tau, s, H) + vmin
            best = v if best is None else min(best, v)
    return best


def generate(case, alpha, beta, N):
    """Exact series of f to R^N: ({(const, m): object array over exponents -rexp..N-rexp}, val_ok)."""
    r = RingQ()
    C = GEN_CASES[case]
    sig, H, rexp = C["sigma"], C["H"], C["rexp"]
    a, b = Fraction(alpha), Fraction(beta)
    kap, lam = (a + b) / 2, (a - b) / 2
    pairs = _pair_terms(C["c"])
    js = sorted({j for j, _ in C["c"]})
    ks = sorted({k for _, k in C["c"]})
    if lam == 0:
        taumax = max(ks) + 2 * H - sig
    else:
        taumax = sig
        while tau_bound(case, taumax + 1) <= N:
            taumax += 1
    Pc, Vc = legendre_PV(max(taumax, 2))
    jmax = max(js)
    nmax = 2 * taumax + 2 * jmax + 4 * H - 2 * sig + 3
    emax = N + 2
    TB = Tables(r, nmax, emax)
    E0 = TB.emin
    RLO, RHI = -rexp, N - rexp
    WR = RHI - RLO + 1
    acc = {k: r.vec(WR) for k in ("1", "g", "l2", "lp")}
    val_ok = True
    kscale = r.vec(TB.W)
    for i in range(TB.W):
        kscale[i] = kap ** (E0 + i)
    for tau in range(sig, taumax + 1):
        bser = {}
        for k in ks:
            v = r.vec(N + 1)
            for m in range(N + 1):
                if lam == 0 and m > 0:
                    break
                mo = eta_moment(k + m, tau, sig, H)
                if mo:
                    v[m] = Fraction(-lam) ** m / factorial(m) * mo
            if np.any(v[:min(_vb(k, tau, sig, H), N + 1)] != 0):
                val_ok = False
            bser[k] = v
        wser = {}
        for (j1, j2), lst in pairs.items():
            w = r.vec(N + 1)
            for k1, k2, cc in lst:
                w = w + r.conv(bser[k1], bser[k2])[:N + 1] * cc
            wser[(j1, j2)] = w
        if all(not np.any(w != 0) for w in wser.values()):
            continue
        Fd, Bd, Cd = tau_polys(tau, sig, H, js, Pc, Vc)
        Stau = {k: r.vec(WR) for k in acc}
        for (j1, j2), w in wser.items():
            if not np.any(w != 0):
                continue
            pcs = g_pieces_gen(r, Fd, Bd, Cd, j1, j2)
            Xh: Dict[str, np.ndarray] = {}
            for basis, poly in pcs.items():
                n = len(poly)
                assert n <= nmax + 1, (n, nmax)
                for comp, tab in TB.T[basis].items():
                    contrib = r.dot(poly, tab[:n])
                    Xh[comp] = contrib if comp not in Xh else Xh[comp] + contrib
            vmin = -(j1 + j2 + 4 * H - 2 * sig + 1)
            for comp, v in Xh.items():
                if np.any(v[: vmin - E0] != 0):
                    val_ok = False
                seg = (v * kscale)[vmin - E0:]
                prod = r.conv(w, seg)
                off = vmin - RLO
                nn = min(len(prod), WR - off)
                if nn > 0:
                    Stau[comp][off:off + nn] = Stau[comp][off:off + nn] + prod[:nn]
        cut = tau_bound(case, tau) - rexp - RLO
        for comp, v in Stau.items():
            if cut > 0 and np.any(v[:min(cut, WR)] != 0):
                val_ok = False
        wt = Fraction(2 * tau + 1) * Fraction(factorial(tau - sig), factorial(tau + sig)) ** 2
        for comp in acc:
            acc[comp] = acc[comp] + Stau[comp] * wt
    em = r.vec(N + 1)
    for n in range(N + 1):
        em[n] = Fraction(-2 * kap) ** n / factorial(n)
    pref = C["pref"](a, b)
    out: Dict[Tuple[str, int], np.ndarray] = {}

    def put(key, v):
        out[key] = out[key] + v if key in out else v
    pl = _prime_logs(kap)
    for comp, v in acc.items():
        full = np.convolve(v, em)[:WR] * pref
        if comp == "lp":                                # ln p = ln kappa + ln R
            put(("1", 1), full)
            for pr, e in pl.items():
                put((f"l{pr}", 0), full * e)
        else:
            put((comp, 0), full)
    return out, val_ok


# Series orders in the fast tier: R^40 at (1, 1) (0.3 s), R^24 at (1, 2) (S 7 s, P 3 s; R^48 costs 51 s for S).
N_FAST = {"S_1_1": 40, "P_1_1": 40, "S_1_2": 24, "P_1_2": 24}
N_SLOW = 64


def _to_dicts(ser) -> Dict[Tuple[str, int], Dict[int, Fraction]]:
    """Array index = power of R (the generator folds the R^rexp prefactor in); keep nonzero entries."""
    comps = {}
    for k, v in ser.items():
        d = {e: Fraction(x) for e, x in enumerate(v) if x != 0}
        if d:
            comps[k] = d
    return comps


@lru_cache(maxsize=None)
def series(key: str, N: int = 0):
    """{(const, m): {exponent: coefficient}} exact to R^N (N = 0: the fast-tier order), val_ok."""
    fam, a, b = CASES_T[key][:3]
    ser, ok = generate(fam, a, b, N or N_FAST[key])
    return _to_dicts(ser), ok


# ---------------------------------------------------------------- operator on log-Laurent series


def _ff_poly(i: int) -> List[int]:
    c = [1]
    for s in range(i):
        new = [0] * (len(c) + 1)
        for e, v in enumerate(c):
            new[e + 1] += v
            new[e] -= s * v
        c = new
    return c


@lru_cache(maxsize=None)
def _ffd_coeffs(i: int, t: int) -> Tuple[int, ...]:
    c = _ff_poly(i)
    for _ in range(t):
        c = [e * c[e] for e in range(1, len(c))] or [0]
    return tuple(c)


def ffd(i: int, t: int, k: int) -> int:
    """t-th k-derivative of the falling factorial k^(i), at integer k."""
    v = 0
    for c in reversed(_ffd_coeffs(i, t)):
        v = v * k + c
    return v


def apply_operator(op, comps, tmax: int) -> Dict[Tuple[str, int, int], Fraction]:
    """Coefficients of L f at R^t (ln R)^m', t <= tmax, per constant.

    R^j D^i [R^k ln^m R] = sum_t C(m, t) (d/dk)^t k^(i) R^{k-i+j} ln^{m-t} R.  Exact for t <= N - r when the
    series is exact to R^N.
    """
    out: Dict[Tuple[str, int, int], Fraction] = {}
    for (c, m), d in comps.items():
        for i, p in enumerate(op):
            for j, cij in enumerate(p):
                if not cij:
                    continue
                for k, fk in d.items():
                    t_ord = k - i + j
                    if t_ord > tmax:
                        continue
                    for t in range(m + 1):
                        w = comb(m, t) * ffd(i, t, k)
                        if w:
                            key = (c, m - t, t_ord)
                            out[key] = out.get(key, Fraction(0)) + w * cij * fk
    return out


def frobenius(op, data: Dict[Tuple[str, int], Dict[int, Fraction]], kmax: int, kmin: int = 0):
    """Frobenius continuation from the values at the R = 0 exponents: (F, roots, residuals at the exponents)."""
    cij = {(i, j): v for i, p in enumerate(op) for j, v in enumerate(p) if v}
    delta = min(j - i for (i, j) in cij)
    k = sp.Symbol("k")
    I = sum(cv * sp.ff(k, i) for (i, j), cv in cij.items() if j - i == delta)
    roots = sorted(int(rt) for rt in sp.roots(sp.Poly(I, k)) if rt.is_integer)
    consts = sorted({c for c, m in data})
    F, resid = {}, {}
    for c in consts:
        mmax = max(m for cc, m in data if cc == c)
        Fc: Dict[int, Dict[int, Fraction]] = {m: {} for m in range(mmax + 1)}
        for kk in range(min(roots + [kmin]), kmax + 1):
            for m in range(mmax, -1, -1):
                lhs, top = Fraction(0), Fraction(0)
                for (i, j), cv in cij.items():
                    k2 = kk + delta + i - j
                    for t in range(mmax - m + 1):
                        w = comb(m + t, t) * ffd(i, t, k2)
                        if not w:
                            continue
                        if k2 == kk and t == 0:
                            top += w * cv
                        else:
                            lhs += w * cv * Fc[m + t].get(k2, Fraction(0))
                if top == 0:
                    Fc[m][kk] = data.get((c, m), {}).get(kk, Fraction(0))
                    resid[(c, m, kk)] = lhs
                else:
                    Fc[m][kk] = -lhs / top
        for m in Fc:
            F[(c, m)] = Fc[m]
    return F, roots, resid


def indicial(op) -> Tuple[bool, sp.Poly]:
    """(Fuchsian at R = 0, indicial polynomial in theta)."""
    th = sp.Symbol("theta")
    r = len(op) - 1
    ords = [next(j for j, v in enumerate(p) if v) for p in op]
    vals = [o - i for i, o in enumerate(ords)]
    v = min(vals)
    ind = sum(sp.Rational(op[i][ords[i]].numerator, op[i][ords[i]].denominator) * sp.ff(th, i)
              for i in range(r + 1) if vals[i] == v)
    return vals[r] == v, sp.Poly(sp.expand(ind), th)


def infinity_rank_and_rates(op):
    """Newton-polygon rank at infinity and the characteristic roots of the top-degree edge (with multiplicity)."""
    r = len(op) - 1
    degs = [len(p) - 1 for p in op]
    w = [degs[j] - j for j in range(r + 1)]
    rank = max(sp.Rational(w[j] - w[r], r - j) for j in range(r))
    D = max(degs)
    lam = sp.Symbol("lam")
    chi = sum(sp.Rational(op[k][D].numerator, op[k][D].denominator) * lam ** k
              for k in range(r + 1) if degs[k] == D)
    P = sp.Poly(chi, lam)
    return rank, dict(sp.roots(P)), P.degree()


def _cval(c, mp):
    return {"1": mp.mpf(1), "g": +mp.euler}.get(c) or mp.log(int(c[1:]))


def frobenius_value(key: str, Rv: str, K: int, dps: int):
    """Value at R = Rv of the frozen operator's Frobenius solution, data = f's values at the R = 0 exponents only."""
    from mpmath import mp
    comps, _ = series(key)
    ex = exponents0(key)
    data = {kk: {k: v for k, v in d.items() if k in ex} for kk, d in comps.items()}
    F, _, resid = frobenius(op_coeffs(key), data, K)
    assert all(v == 0 for v in resid.values())
    with mp.workdps(dps):
        x = mp.mpf(Rv)
        tot = mp.mpf(0)
        for (c, m), d in F.items():
            tot += _cval(c, mp) * mp.log(x) ** m * mp.fsum(mp.mpf(v.numerator) / v.denominator * x ** k
                                                         for k, v in d.items() if v)
        last = max(abs(mp.mpf(v.numerator) / v.denominator * x ** k)
                   for (c, m), d in F.items() for k, v in d.items() if k >= K - 3 and v)
        return +tot, +last


# hydrogenic orbitals for the production routes, as ep_gen.REF_ORBITALS: (ZA, oa, ob, ZB, oc, od); the hydrogenic
# 2p rate is Z_A / 2, so Z_A = 2a.
REF_ORBITALS = {
    "2psig": lambda a, b: (2 * a, (2, 1, 0), (1, 0, 0), b, (2, 1, 0), (1, 0, 0)),
    "2ppi": lambda a, b: (2 * a, (2, 1, 1), (1, 0, 0), b, (2, 1, -1), (1, 0, 0)),
}


@lru_cache(maxsize=None)
def closed_form_eq(key: str):
    """f(R) at equal exponents from geovac.qfd_core.exchange_closed_form, symbolic in R (the tau sum terminates;
    the terms past tau = 3 are checked to vanish identically)."""
    from geovac import qfd_core as Q
    fam, a, b = CASES_T[key][:3]
    assert a == b
    Rp = sp.Symbol("Rp", positive=True)
    orb = REF_ORBITALS[fam](sp.Rational(a.numerator, a.denominator), sp.Rational(b.numerator, b.denominator))
    terms = Q.exchange_closed_form(*orb, Rp, tau_max=6, return_terms=True)
    assert all(t == 0 for t in terms[-2:])
    return Rp, sp.Add(*terms)


def ode_residual_closed_form(key: str, Rv: str, dps: int = 60):
    """(|L f|, sum_i |p_i f^(i)|) at R = Rv, f the production closed form; derivatives by mpmath.diff."""
    from mpmath import mp
    Rp, expr = closed_form_eq(key)
    op = op_coeffs(key)
    with mp.workdps(dps):
        fn = sp.lambdify(Rp, expr, "mpmath")
        x = mp.mpf(Rv)
        ders = list(mp.diffs(fn, x, len(op) - 1))
        terms = []
        for i, p in enumerate(op):
            pv = mp.fsum(mp.mpf(c.numerator) / c.denominator * x ** j for j, c in enumerate(p) if c)
            terms.append(pv * ders[i])
        return abs(mp.fsum(terms)), mp.fsum(abs(t) for t in terms)


# ================================================================ tests


@pytest.mark.parametrize("key", KEYS)
def test_operator_shape(key):
    """Order 4 / 4 / 9 / 9, degree 16 / 9 / 34 / 25; every p_i nonzero.

    Rejects: the 1s operator substituted (order 4 degree 8 at (1, 1), order 9 degree 16 at (1, 2)), the other
    family's operator, a truncated or padded coefficient list.
    """
    op = op_coeffs(key)
    order, degree = CASES_T[key][4], CASES_T[key][5]
    assert (len(op) - 1, max(len(p) - 1 for p in op)) == (order, degree)
    assert {k: (CASES_T[k][4], CASES_T[k][5]) for k in KEYS} == {
        "S_1_1": (4, 16), "P_1_1": (4, 9), "S_1_2": (9, 34), "P_1_2": (9, 25)}
    assert all(any(p) for p in op)


@pytest.mark.parametrize("key", KEYS)
def test_operator_annihilates_exact_series(key):
    """L f = 0 exactly over Q, every (constant, ln-power) component, every order t <= N - r (N = 40 at (1, 1),
    24 at (1, 2); the slow tier repeats at N = 64).

    Rejects: any perturbed operator coefficient c_ij with j - i <= N - r (it enters at t = j - i, from f's R^0 term;
    at (1, 1) that is every coefficient, at (1, 2) those with j - i <= 15 -- the rest in the slow tier), the other
    family's or the other pair's operator, and a series generator with a wrong sigma branch, Neumann weight, density
    polynomial or exponent pair.  Non-vacuous: the series has the components the record names (incl. ln 3 at
    kappa = 3/2) and the residual system has >= 3 x 15 checked entries.
    """
    comps, ok = series(key)
    assert ok, "per-tau valuation check failed"
    a, b = CASES_T[key][1:3]
    kap = (a + b) / 2
    primes = {(f"l{q}", 0) for q in _prime_logs(kap) if q != 2}
    assert set(comps) == {("1", 0), ("1", 1), ("g", 0), ("l2", 0)} | primes
    op = op_coeffs(key)
    tmax = N_FAST[key] - (len(op) - 1)
    res = apply_operator(op, comps, tmax)
    assert len(res) >= 3 * 15
    bad = {k: v for k, v in res.items() if v != 0}
    assert not bad, sorted(bad)[:5]


def test_operators_are_case_specific():
    """Negative control: each operator FAILS on the other family's series at the same pair and on the same family's
    series at the other pair.  Rejects an annihilation check that passes anything."""
    for k1, k2 in (("S_1_1", "P_1_1"), ("P_1_1", "S_1_1"), ("S_1_1", "S_1_2"), ("P_1_2", "P_1_1"),
                   ("S_1_2", "P_1_2"), ("P_1_2", "S_1_2")):
        s, _ = series(k2)
        assert any(v != 0 for v in apply_operator(op_coeffs(k1), s, 14).values()), (k1, k2)


@pytest.mark.parametrize("key", KEYS)
def test_R0_fuchsian_exponents(key):
    """R = 0 regular singular; exponents {-3, 0, 4, 11} at (1, 1), {-3 (x2), -1, 0, 1, 2, 4, 6, 11} at (1, 2).

    Rejects: the 1s exponents (top exponent 9, lowest -1), a simple -3 at (1, 2), any other set, a non-Fuchsian point.
    """
    fuchs, ind = indicial(op_coeffs(key))
    assert fuchs
    assert ind.degree() == len(op_coeffs(key)) - 1
    assert {int(k): int(m) for k, m in sp.roots(ind).items()} == exponents0(key)


@pytest.mark.parametrize("key", KEYS)
def test_R0_frobenius_reproduces_series(key):
    """Frobenius from f's values at the exponents reproduces the generator's series to R^N in every component,
    and every equation at an exponent is satisfied.  Rejects an operator whose recursion disagrees with the series
    between the data orders, and one for which f's data are inconsistent."""
    comps, _ = series(key)
    ex = exponents0(key)
    data = {kk: {k: v for k, v in d.items() if k in ex} for kk, d in comps.items()}
    kmax = N_FAST[key]
    F, roots, resid = frobenius(op_coeffs(key), data, kmax)
    assert roots == sorted(ex)
    assert all(v == 0 for v in resid.values()), {k: v for k, v in resid.items() if v}
    assert min(min(d) for d in comps.values()) == 0          # f is finite at R = 0; exponent -3 carries no data
    for kk, d in comps.items():
        for k in range(-3, kmax + 1):
            assert F[kk].get(k, 0) == d.get(k, 0), (kk, k)


@pytest.mark.parametrize("key", KEYS)
def test_lnR_onset_and_coefficient_from_series(key):
    """f's ln R part starts at R^11 ln R with coefficient 704, 16, 22528, 512 (/1157625) for S11, P11, S12, P12.

    Rejects: the 1s onset R^9, any other coefficient (e.g. the 1s (8/3375) a^5 b^5), an extra ln R component.
    """
    comps, _ = series(key)
    assert {m for _, m in comps} == {0, 1}
    assert [c for c, m in comps if m == 1] == ["1"]
    ln = comps[("1", 1)]
    assert min(ln) == 11
    assert ln[11] == LN_COEF[key]
    assert LN_COEF[key] == {"S_1_1": Fraction(704, 1157625), "P_1_1": Fraction(16, 1157625),
                            "S_1_2": Fraction(22528, 1157625), "P_1_2": Fraction(512, 1157625)}[key]


def test_lnR_coefficient_ratio_32():
    """(1,2)/(1,1) = 32 = 2^5 in each family: the R^11 ln R coefficient scales as b^5 at a = 1.

    Rejects: ratio 2^7 (an a^5 b^7-type scaling), 2^3 (the prefactor's b^3 alone), and any family mismatch."""
    assert LN_COEF["S_1_2"] / LN_COEF["S_1_1"] == 32
    assert LN_COEF["P_1_2"] / LN_COEF["P_1_1"] == 32
    s11, _ = series("S_1_1")
    p11, _ = series("P_1_1")
    s12, _ = series("S_1_2")
    p12, _ = series("P_1_2")
    assert s12[("1", 1)][11] / s11[("1", 1)][11] == 32
    assert p12[("1", 1)][11] / p11[("1", 1)][11] == 32


@pytest.mark.parametrize("key", KEYS)
def test_lnR_coefficient_forced_by_operator(key):
    """The R^11 ln R coefficient is FORCED by the operator from f's non-log data below R^11.

    At the exponent 11 the ln-free equation is affine in the R^11 ln R coefficient x; solving it gives the record's
    value.  Independent of the generator's own ln R output; non-zero x means the local solution space at R = 0 is
    not log-free.  Rejects an operator that carries a different (or no) log coefficient, and a wrong expected value.
    """
    comps, _ = series(key)
    ex = exponents0(key)
    base = {kk: {k: v for k, v in d.items() if k in ex and k < 11} for kk, d in comps.items() if kk[1] == 0}
    vals = []
    for x in (Fraction(0), Fraction(1)):
        d = dict(base)
        d[("1", 1)] = {11: x}
        _, _, resid = frobenius(op_coeffs(key), d, 11)
        vals.append(resid[("1", 0, 11)])
    assert vals[1] != vals[0]
    x = -vals[0] / (vals[1] - vals[0])
    assert x == LN_COEF[key]
    assert x != 0


@pytest.mark.parametrize("key", KEYS)
def test_infinity_rank_one_rates(key):
    """Irregular rank exactly 1 at infinity; every behaviour on the top-degree edge, rates {0, 2, -2 (x2)} at (1, 1)
    and {-3, 4, 2, 0 (x2), -2 (x2), -4 (x2)} at (1, 2); no a - b rate.

    Rejects: an inserted +-(a - b) = +-1 rate, unsplit rates at (1, 2), a wrong multiplicity, a rank other than 1.
    """
    a, b = CASES_T[key][1:3]
    rank, roots, deg = infinity_rank_and_rates(op_coeffs(key))
    assert rank == 1
    assert deg == len(op_coeffs(key)) - 1
    assert roots == rates_inf(a, b)
    assert sp.Integer(1) not in roots and sp.Integer(-1) not in roots


@pytest.mark.parametrize("key", KEYS_EQ)
@pytest.mark.parametrize("Rv", ["1.4", "3.015"])
def test_ode_residual_on_production_closed_form(key, Rv):
    """INDEPENDENT (no generator, no series): L f = 0 at R = Rv with f = geovac.qfd_core.exchange_closed_form, the
    exact symbolic production route (terminating tau sum at equal exponents); relative residual < 1e-40 at dps 60.

    Rejects: an operator for the other family, an operator off in any coefficient by a relative amount above
    ~1e-35, and a closed form with the wrong orbitals (e.g. 2p-sigma for 2p-pi).
    """
    res, scale = ode_residual_closed_form(key, Rv)
    assert scale > 0
    assert res / scale < 1e-40, (res, scale)


@pytest.mark.parametrize("key", KEYS_EQ)
def test_frobenius_matches_production_closed_form(key):
    """INDEPENDENT check of the series route: the frozen operator's Frobenius solution from the generator's data at
    the four exponents equals the production closed form at R = 1.4 to 1e-40.

    Rejects: a series generator with a wrong prefactor (e.g. /32 for /64), sigma branch or orbital family, and
    Frobenius data at the wrong exponents."""
    from mpmath import mp
    vf, last = frobenius_value(key, "1.4", K=150, dps=60)
    assert last < mp.mpf(10) ** -50
    Rp, expr = closed_form_eq(key)
    with mp.workdps(60):
        vc = sp.lambdify(Rp, expr, "mpmath")(mp.mpf("1.4"))
        assert abs(vf - vc) < mp.mpf(10) ** -40, (vf, vc)


def test_frobenius_matches_production_exchange_hp_fast():
    """INDEPENDENT at (1, 2), fast: P_1_2's Frobenius solution (data at the nine exponents only) vs the live
    production route geovac.qfd_core.exchange_hp at R = 0.7, tau <= 7, dps 30 (last tau term 2.3e-15, checked
    < 1e-14), to 1e-13.  S_1_2 and R = 1.4 are in the slow tier (exchange_hp costs 46-190 s there).

    Rejects: an operator or generator data whose continuation is not the exchange integral (the S_1_2 operator on
    P data, a prefactor off by a factor, the (1, 1) data) -- all move the value at >= 1e-3 relative."""
    from mpmath import mp
    from geovac import qfd_core as Q
    fam, a, b = CASES_T["P_1_2"][:3]
    vf, last = frobenius_value("P_1_2", "0.7", K=100, dps=50)
    assert last < mp.mpf(10) ** -40
    orb = REF_ORBITALS[fam](sp.Rational(a.numerator, a.denominator), sp.Rational(b.numerator, b.denominator))
    vh, per = Q.exchange_hp(*orb, sp.Rational("0.7"), tau_max=7, dps=30)
    assert abs(per[-1]) < mp.mpf(10) ** -14
    assert abs(vf - vh) < mp.mpf(10) ** -13, (vf, vh)


# ---------------------------------------------------------------- slow tier


@pytest.mark.slow
@pytest.mark.parametrize("key", KEYS)
def test_operator_annihilates_exact_series_slow(key):
    """The exact check over Q to R^64: orders up to 64 - r, past every coefficient's entry order j - i <= 34, so
    a perturbation of ANY coefficient fires (the fast tier at (1, 2) reaches j - i <= 15 only)."""
    comps, ok = series(key, N_SLOW)
    assert ok
    op = op_coeffs(key)
    tmax = N_SLOW - (len(op) - 1)
    assert tmax >= max(len(p) - 1 for p in op)
    res = apply_operator(op, comps, tmax)
    assert len(res) >= 3 * 50
    bad = {k: v for k, v in res.items() if v != 0}
    assert not bad, sorted(bad)[:5]


@pytest.mark.slow
@pytest.mark.parametrize("key", KEYS_UNEQ)
def test_frobenius_matches_production_exchange_hp(key):
    """INDEPENDENT at (1, 2): the frozen operator's Frobenius solution from f's data at the nine exponents vs the
    production route geovac.qfd_core.exchange_hp (per-tau closed forms of the xi and eta integrals, mpmath sum) at
    R = 1.4, to 1e-20 (the omitted tau tail is checked below 1e-24).

    Rejects: an operator or data whose Frobenius continuation is not the exchange integral, including the other
    family's operator and a generator prefactor off by a factor."""
    from mpmath import mp
    from geovac import qfd_core as Q
    fam, a, b = CASES_T[key][:3]
    vf, last = frobenius_value(key, "1.4", K=150, dps=60)
    assert last < mp.mpf(10) ** -40
    orb = REF_ORBITALS[fam](sp.Rational(a.numerator, a.denominator), sp.Rational(b.numerator, b.denominator))
    vh, per = Q.exchange_hp(*orb, sp.Rational("1.4"), tau_max=TAU_HP, dps=30)
    assert abs(per[-1]) < mp.mpf(10) ** -24
    assert abs(vf - vh) < mp.mpf(10) ** -20, (vf, vh)


TAU_HP = 12


# TIMINGS (2026-11-17, this machine): fast tier 37 tests in ~27 s (series: S12 to R^24 7 s, P12 3 s, (1, 1) < 0.5 s;
# live exchange_hp P12 R = 0.7 tau <= 7: 11 s; closed-form ODE residuals 1-3 s each).  Slow tier 6 tests, ~285 s
# (--slow): exact annihilation to R^64, S12 88 s, P12 28 s; live exchange_hp R = 1.4 tau <= 12, S12 125 s, P12 41 s.
# Measured here: closed-form ODE relative residuals 2e-61 .. 8e-60 (bar 1e-40); Frobenius vs closed form at R = 1.4,
# 7e-61 and 2e-60 (bar 1e-40); Frobenius vs exchange_hp agrees to the reference's own tau tail.
#
# FIRE-TEST PLANTS (python debug/qa/fire_test.py tests/test_paper58_exchange_l1.py --plant "OLD=>NEW" -k <sel>),
# 2026-11-17, a separate step after the guard; 37 runs, 35 fired; the 2 that did not are by design (below).
# ANN = "annihilates_exact_series and not slow".
#     S11 p_0 coeff +1                    "2560*R**14" -> "2561*R**14"                     -k "ANN and S_1_1"            FIRED
#     same, independent route             same                                             -k "ode_residual and S_1_1"   FIRED
#     P11 p_0 coeff +1                    "96*R**7 + 704" -> "... 705"                     -k "ANN and P_1_1"            FIRED
#     same, independent route             same                                             -k "ode_residual and P_1_1"   FIRED
#     P11 p_4 leading coeff +1            "R**3*(96*R**6" -> "R**3*(97*R**6"               -k "ode_residual and P_1_1"   FIRED
#     S12 p_0 constant +1                 "...370384000)" -> "...370384001)"               -k "ANN and S_1_2"            FIRED
#     P12 p_0 constant +1                 "...788104000)" -> "...788104001)"               -k "ANN and P_1_2"            FIRED
#     same, independent route             same                                             -k exchange_hp_fast           FIRED
#     S11 given the P11 operator          CASES_T OP_S_1_1 -> OP_P_1_1                     -k "S_1_1 and (ANN or ode_residual)"  FIRED
#     same                                same                                             -k "operator_shape and S_1_1" FIRED
#     S12 given the S11 operator          CASES_T OP_S_1_2 -> OP_S_1_1                     -k "S_1_2 and (shape or ANN)" FIRED
#     P12 given the S12 operator          CASES_T OP_P_1_2 -> OP_S_1_2                     -k exchange_hp_fast           FIRED
#     exponent 11 -> 9 at (1, 1)          EXPONENTS_EQ "11: 1" -> "9: 1"                   -k "fuchsian and 1_1"         FIRED
#     -3 simple at (1, 2)                 "{-3: 2," -> "{-3: 1,"                           -k "fuchsian and 1_2"         FIRED
#     exponent 6 -> 5 (Frobenius data)    "6: 1, 11: 1}" -> "5: 1, 11: 1}"                 -k "frobenius_reproduces and 1_2"  FIRED
#     ln coeff 704 -> 705                 LN_COEF S_1_1                                    -k "lnR_onset and S_1_1"      FIRED
#     same, operator-forced route         same                                             -k "lnR_coefficient_forced and S_1_1"  FIRED
#     ln coeff 512 -> 1024 (ratio 64)     LN_COEF P_1_2                                    -k ratio_32                   FIRED
#     same, operator-forced route         same                                             -k "lnR_coefficient_forced and P_1_2"  FIRED
#     rate a - b inserted                 "q(-(a + b)): 1" -> "q(a - b): 1"                -k "infinity and 1_2"         FIRED
#     (1, 1) rate multiplicities swapped  "q(2 * a): 1, q(-2 * a): 2" -> "2 ... 1"         -k "infinity and 1_1"         FIRED
#     Neumann weight 2 tau + 1 -> + 3     generator                                        -k ANN                        FIRED
#     2psig density (1, 3) sign flipped   GEN_CASES "(1, 3): -1" -> "(1, 3): 1"            -k "S_1_1 and (ANN or frobenius_matches_production_closed_form)"  FIRED
#     2ppi prefactor /64 -> /32           GEN_CASES                                        -k "frobenius_matches_production_closed_form and P_1_1"  FIRED
#     same                                same                                             -k "lnR_onset and P"          FIRED
#     same                                same                                             -k exchange_hp_fast           FIRED
#     kappa, lambda swapped in generator  "(a + b) / 2, (a - b) / 2" swapped               -k "ANN and 1_2"              FIRED
#     ln-R derivative of k^(i) dropped    "[e * c[e] ..." -> "[c[e] ..."                   -k ANN                        FIRED
#     closed form: 2p-sigma orbitals for P  REF_ORBITALS 2ppi -> (2,1,0),(2,1,0)           -k "ode_residual and P_1_1"   FIRED
#     closed form: Z_A = a (2p rate a/2)  REF_ORBITALS 2psig "2 * a" -> "1 * a"            -k "ode_residual and S_1_1"   FIRED
#     exchange_hp tau_max 7 -> 3          test_frobenius_matches_production_exchange_hp_fast  -k exchange_hp_fast        FIRED
#     apply_operator output zeroed        "+ w * cij * fk" -> "+ 0 * w * cij * fk"         -k case_specific              FIRED
#     S12 c_{0,32} +1, slow tier          "445510172085608448000*R**32" -> "...001"        -k "annihilates_exact_series_slow and S_1_2"  FIRED
#     P12 c_{0,23} +1, slow tier          "384*(1325049062400*R**23" -> "...401"           -k "annihilates_exact_series_slow and P_1_2"  FIRED
#     2psig prefactor /32 -> /16, slow    GEN_CASES                                        -k "exchange_hp and S_1_2"    FIRED
#   Did not fire, by design:
#     apply_operator output zeroed, -k ANN: a harness that returns all-zero residuals passes any annihilation check;
#       that is what test_operators_are_case_specific is for, and it fires on the same plant (above).
#     S12 c_{0,32} +1, fast tier: the coefficient enters at order 32 > 24 - 9 = 15, past the fast series; the slow
#       tier (R^64) fires on it (above).  Fast-tier coverage at (1, 2) is coefficients with j - i <= 15.
