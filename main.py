from tkinter import messagebox
import matplotlib.pyplot as plt
import re
import tkinter as Tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

class gff_bestand():
    def __init__(self,bestandsnaam):
        self.__bestandsnaam = bestandsnaam
        self.bestand_in_lijst()
        self.lijsten_aanmaken()


    def bestand_in_lijst(self):
        bestand = open(self.__bestandsnaam)
        self.__inhoud_gff = []
        for regel in bestand:
            regel = regel.replace("\n", "")  # Haalt de enter eruit
            lijst = regel.split("\t")  # splits de lijst op elke tab
            self.__inhoud_gff.append(lijst)
        bestand.close()

    def lijsten_aanmaken(self):
        self.__cds = []
        self.__gene = []
        self.__trna = []
        self.__exon = []
        self.__anders = []
        for lijst in self.__inhoud_gff:
            if lijst[2] == "gene":
                self.__gene.append(lijst)
            elif lijst[2] == "CDS":
                self.__cds.append(lijst)
            elif lijst[2] == "tRNA":
                self.__trna.append(lijst)
            elif lijst[2] == "exon":
                self.__exon.append(lijst)
            else:
                self.__anders.append(lijst)

    def getlijsten(self):
        return self.__cds, self.__gene, self.__trna, self.__exon, \
               self.__anders, self.__inhoud_gff

class gffb_bestand():
    def __init__(self,bestandsnaam):
        self.__bestandsnaam = bestandsnaam
        self.bestand_opdelen()
        self.informatie_halen_uit_bestand2()
        self.dictonaries_aanmaken()

    
    def bestand_opdelen(self):
        self.__bestand1_tekst = open("bestand1", "w")
        self.__bestand2_informatie = open("gene+cds.txt", "w")
        self.__bestand3_sequentie = open("sequentie.txt", "w")

        bestand = open(self.__bestandsnaam, "r")
        regel_nummer = 0
        for regel in bestand:
            if regel_nummer >= 0 and regel_nummer <= 45:
                self.__bestand1_tekst.write(regel)
            elif regel_nummer <= 49099 and regel_nummer >= 46:
                self.__bestand2_informatie.write(regel)
            elif regel_nummer >= 49100:
                self.__bestand3_sequentie.write(regel)
            regel_nummer += 1
        self.__bestand1_tekst.close()
        self.__bestand2_informatie.close()
        self.__bestand3_sequentie.close()

    def informatie_halen_uit_bestand2(self):
        """

        :return: lijst_met_gene_cds - 2d lijst met inhoudt van bestand gene+cds.txt
        - [["gene 517..1878", /gene="dnaA", ...] [CDS 517..1878, /gene="dnaA",...],
          ["gene ..., n]]
        """
        bestand = open(self.__bestand2_informatie.name)
        # bestand = open("test.txt")
        self.__lijst_met_gegevens = []
        per_gegeven = []
        per_kopje_in_bestand = ""
        aantal = 0
        for regel in bestand:
            while regel.startswith(" "):  # Haalt alle spaties weg dat bij het
                regel = regel[1:]  # begin van de regel staat
            regel = regel.replace("\n", "")
            # zorgt ervoor dat elke keer wanneer "gene" of "CDS"staat er een nieuw
            # lijst wordt aangemaakt en de oude wordt toegevoegt aan de 2d lijst
            if regel.startswith("gene") or regel.startswith("CDS") or regel.startswith("rRNA   "):
                per_gegeven.append(per_kopje_in_bestand)
                self.__lijst_met_gegevens.append(per_gegeven)
                per_gegeven = []
                per_gegeven.append(regel)
                per_kopje_in_bestand = ""
            else:
                # Zorgt ervoor dat kopjes die bestaan uit meerdere regels
                # samen als een regel worden toegevoegt in een lijst
                if regel.startswith("/"):
                    if per_kopje_in_bestand != "":
                        per_gegeven.append(per_kopje_in_bestand)
                        per_kopje_in_bestand = ""
                    per_kopje_in_bestand += regel
                else:
                    per_kopje_in_bestand += regel
        per_gegeven.append(per_kopje_in_bestand)
        self.__lijst_met_gegevens.append(per_gegeven)
        #    aantal += 1
        #     if aantal == 50:
        #         break

        del self.__lijst_met_gegevens[0]

    def dictonaries_aanmaken(self):
        self.__gene = {}
        self.__cds = {}
        self.__rrna = {}
        id_nummer_gene = 1
        id_nummer_cds = 1
        id_nummer_rrna = 1
        for lijst in self.__lijst_met_gegevens:
            if lijst[0].startswith("gene  "):
                self.__gene[id_nummer_gene] = lijst
                id_nummer_gene += 1
            elif lijst[0].startswith("CDS  "):
                self.__cds[id_nummer_cds] = lijst
                id_nummer_cds += 1
            else:
                self.__rrna[id_nummer_rrna] = lijst
                id_nummer_rrna += 1

    def getdiconaries(self):
        return self.__gene, self.__cds, self.__rrna

class Grafieken():

    def __init__(self):
        self.__grafieksoort = input("Welk grafiek wil je zien?\n\n"
                                    "Kies uit:\n- aantal sequenties"
                                    "+/- strand\n - ")


    def labels(x,y):
        for i in range(len(x)):
            plt.text(i, y[i], y[i])


class Taartdiagram(Grafieken):

    def __init__(self, gegevens):
        self.__gegevens = gegevens
        #self.setgegevens_strands()


    def setgegevens_strands(self):
        self.__title = "Aantal sequenties op +/- strand"
        self.__labels = "+","-"
        self.__grote = [0,0]
        for lijst in self.__gegevens:
            if lijst[6] == "+":
                self.__grote[0] += 1
            elif lijst[6] == "-":
                self.__grote[1] += 1

    def setgegevens_gc(self):
        self.__title = "GC/AT %"
        self.__labels = "GC%", "AT%"
        gc = round((self.__gegevens.count("c") + self.__gegevens.count("g")) / len(self.__gegevens),2) *100
        at = round((self.__gegevens.count("a") + self.__gegevens.count("t")) / len(self.__gegevens),2) *100
        self.__grote = [gc,at]

    def maakgrafiek(self):
        figuur, vormgeving = plt.subplots()
        vormgeving.pie(self.__grote, labels=self.__labels)
        plt.title(self.__title)
        plt.text(0,0.6, self.__grote[0])
        plt.text(0,-0.6, self.__grote[1])
        plt.show()


class Staafdiagram(Grafieken):

    def __init__(self, gegevens):
        self.__gegevens = gegevens

    def labels(self):
        for i in range(len(self.__x)):
            plt.text(i, self.__y[i], self.__y[i])

    def setgegevenstypes(self):
        self.__labels = {}
        for lijst in self.__gegevens:
            if lijst[2] not in self.__labels:
                self.__labels[lijst[2]] = 1
            else:
                self.__labels[lijst[2]] += 1
        self.__x = []
        self.__y = []
        for key in self.__labels.keys():
            self.__x.append(key)
        for value in self.__labels.values():
            self.__y.append(value)
        self.labels()


    def setgegevenshypothetical_proteins(self):
        aantal_wel = 0
        aantal_niet = 0
        for k,v in self.__gegevens.items():
            for item in v:
                if item == "/note=\"conserved hypothetical protein\"":
                    aantal_wel += 1
                else:
                    aantal_niet += 1
        self.__x = ["aantal hypothetical protein","anders"]
        self.__y = [aantal_wel,aantal_niet]
        self.labels()

    def maakgrafiek(self):
        plt.bar(self.__x,self.__y)
        plt.title("grafiek met hoevaak elke type voorkomt")
        plt.xlabel("Type")
        plt.ylabel("aantal")
        plt.grid()
        plt.show()


class gui_maken():
    def __init__(self,patronen,gff_totale_inhoud,gffb_cds):
        self.__patronen = patronen
        self.__gff_totale_inhoud = gff_totale_inhoud
        self.__gffb_cds = gffb_cds
        self.setgegevenstypes()
        self.setgegevenshypothetical_proteins()
        self.settabel()
        self.setgegevens_strands()

        self.gui_startscherm()


    def setgegevenstypes(self):
        self.__labels = {}
        for lijst in self.__gff_totale_inhoud:
            if lijst[2] not in self.__labels:
                self.__labels[lijst[2]] = 1
            else:
                self.__labels[lijst[2]] += 1
        self.__x = []
        self.__y = []
        for key in self.__labels.keys():
            self.__x.append(key)
        for value in self.__labels.values():
            self.__y.append(value)

    def settabel(self):
        self.__tabel = []
        self.__tabel.append(["nummer","eiwit","protein_id","product","+/- strand"])
        for patroon in self.__patronen:
            self.__tabel.append(patroon)


    def gui_startscherm(self):
        startscherm = Tk.Tk()
        self.__knop = Tk.Button(heigh=2, master=startscherm, text="klik hier voor tabel patronen",command=self.tabel_patronen)
        self.__knop2 = Tk.Button(heigh=2, master=startscherm, text="klik hier voor grafiek met verschillende type ",command=self.grafiek_type)
        self.__knop3 = Tk.Button(heigh=2, master=startscherm,text="klik hier voor grafiek met hypotheticale proteins ",command=self.grafiek_hypothetical_proteins)
        self.__knop4 = Tk.Button(heigh=2, master=startscherm,text="klik hier voor grafiek met % +/- strands ",command=self.maak_grafiek_strands)
        self.__knop.pack()
        self.__knop2.pack()
        self.__knop3.pack()
        self.__knop4.pack()
        Tk.mainloop()

    def tabel_patronen(self):
        tabelscherm =Tk.Toplevel()
        tabelscherm.title("gevonden patronen")
        aantal_rijen = len(self.__tabel)
        aantal_kolommen = len(self.__tabel[0])
        for rij in range(aantal_rijen):
            for kolom in range(aantal_kolommen):
                tabel =Tk.Entry(tabelscherm,fg="black",width=50)
                tabel.grid(row=rij,column=kolom)
                tabel.insert(Tk.END, self.__tabel[rij][kolom])
        knop = Tk.Button(text="return", master=tabelscherm, command=tabelscherm.destroy)
        knop.grid()

    def grafiek_type(self):
        grafiek_type =Tk.Toplevel()
        keys = []
        data = []
        for k, v in self.__labels.items():
            data.append(v)
            keys.append(k)
        f = Figure(figsize=(10, 8), dpi=100)
        ax = f.add_subplot(111)
        width = .8
        rects1 = ax.bar(keys, data, width)
        figuur = FigureCanvasTkAgg(f, master=grafiek_type)
        figuur.draw()
        figuur.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        quit_button = Tk.Button(heigh=2, width=10, text="return",
                    command=grafiek_type.destroy, master=grafiek_type)
        quit_button.pack()

    def setgegevenshypothetical_proteins(self):
        aantal_wel = 0
        aantal_niet = 0
        for k,v in self.__gffb_cds.items():
            for item in v:
                if item == "/note=\"conserved hypothetical protein\"":
                    aantal_wel += 1
                else:
                    aantal_niet += 1
        self.__x = ["aantal hypothetical protein","anders"]
        self.__y = [aantal_wel,aantal_niet]
        #self.labels()

    def grafiek_hypothetical_proteins(self):
        grafiek_hypothetical_proteins = Tk.Toplevel()
        f = Figure(figsize=(10, 8), dpi=100)
        ax = f.add_subplot(111)
        width = .8
        rects1 = ax.bar(self.__x, self.__y, width)
        figuur = FigureCanvasTkAgg(f, master=grafiek_hypothetical_proteins)
        figuur.draw()
        figuur.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        quit_button = Tk.Button(heigh=2, width=10, text="return",
                                command=grafiek_hypothetical_proteins.destroy,
                                master=grafiek_hypothetical_proteins)
        quit_button.pack()

    def setgegevens_strands(self):
        self.__grote = [0, 0]
        for lijst in self.__gff_totale_inhoud:
            if lijst[6] == "+":
                self.__grote[0] += 1
            elif lijst[6] == "-":
                self.__grote[1] += 1
        totaal = self.__grote[0] + self.__grote[1]
        # zodat het in het rondje past
        self.__grote[0] = (self.__grote[0]/totaal) *1000
        self.__grote[1] = (self.__grote[1] / totaal) * 1000

    def maak_grafiek_strands(self):
        grafiek_strands = Tk.Toplevel()
        upper_frame = Tk.Frame(grafiek_strands)
        bottom_frame = Tk.Frame(grafiek_strands)
        upper_frame.pack()
        bottom_frame.pack()
        label = Tk.Label(upper_frame, text='percentage +/-strands')
        legenda = Tk.Label(bottom_frame, text=("legenda:\n -red: +\n-blue: -"))
        label.pack()
        legenda.pack()
        pie = Tk.Canvas(upper_frame, width=154, height=154)
        pie.pack()
        pie.create_arc((2, 2, 152, 152), fill="red", start=0, extent=(360*self.__grote[0])/1000)
        pie.create_arc((2, 2, 152, 152), fill="blue",start= (360*self.__grote[0])/1000, extent=(360*self.__grote[1])/1000)
        quit_button = Tk.Button(heigh=2, width=10, text="return",
                                command=grafiek_strands.destroy,
                                master=bottom_frame)
        quit_button.pack()



class Patronen():
    def __init__(self,cds,gffb_gene,gffb_cds):
        self.__gff_cds = cds
        self.__gffb_cds = gffb_cds
        self.__gffb_gene = gffb_gene
        self.vindpatronen()
        self.vindprotein_id()
        self.vindproduct()
        self.vindstrand()



    def vindpatronen(self):
        self.__gevonden_patronen = []
        aantal = 1
        for key, value in self.__gffb_cds.items():
            patroon = []
            for i in value:
                if i.startswith("/translation"):
                    x = bool(re.match(
                        ".*[ST]G[LIVMFYW]{3}[GN][A-Z]{2}T[LIVM][A-Z]T[A-Z]{2}H.*",
                        i))
                    y = bool(
                        re.match(".*T[A-Z]{2}[GC][NQ]SGS[A-Z][LIVM][FY].*", i))
                    if x == True or y == True:
                        patroon.append(aantal)
                        i = i.replace("/", "")
                        patroon.append(i)
            if patroon != []:
                self.__gevonden_patronen.append(patroon)
                aantal += 1

    def vindprotein_id(self):
        for patroon in self.__gevonden_patronen:
            value = self.__gffb_cds[patroon[0]]
            for i in value:
                if i.startswith("/protein_id="):
                    i = i.replace("/", "")
                    patroon.append(i)

    def vindproduct(self):
        for patroon in self.__gevonden_patronen:
            value = self.__gffb_cds[patroon[0]]
            for i in value:
                if i.startswith("/product="):
                    i = i.replace("/","")
                    patroon.append(i)

    def vindstrand(self):
        for patroon in self.__gevonden_patronen:
            key = patroon[0]
            value = self.__gff_cds[key]
            strand = value[6]
            patroon.append(strand)

    def getpatronen(self):
        return self.__gevonden_patronen


def grafieken(inhoud_gff,sequentie_genoom,cds):
    grafieksoort = input("Welk grafiek wil je zien?\n\nKies uit:\n"
                         "- aantal sequenties +/- strand\n- types\n"
                         "- GC%\n- aantal hypothetical proteins\n")
    if grafieksoort == "aantal sequenties +/- strand":
        g = Taartdiagram(inhoud_gff)
        g.setgegevens_strands()
        print (g.maakgrafiek())
    elif grafieksoort == "types":
        g = Staafdiagram(inhoud_gff)
        g.setgegevenstypes()
        print (g.maakgrafiek())
    elif grafieksoort == "GC%":
        g = Taartdiagram(sequentie_genoom)
        g.setgegevens_gc()
        print (g.maakgrafiek())
    elif grafieksoort == "aantal hypothetical proteins":
        g = Staafdiagram(cds)
        g.setgegevenshypothetical_proteins()
        print (g.maakgrafiek())

    else:
        print("Verzoek niet gevonden! probeer het opnieuw")
        grafieken(inhoud_gff)

def main():
    bestandsnaam_gff = "gff_bestand.txt"
    bestandsnaam_gbff = "gbff.txt"

    #inhoud_gff = informatieinlezen_gff(bestandsnaam_gff)
    #gene_cds_bestand, sequentie_bestand = gbff_bestand_opdelen(bestandsnaam_gbff)

    gff = gff_bestand(bestandsnaam_gff)
    gff_cds, gff_gene, gff_trna, ggf_exon, gff_anders, gff_totale_inhoud = gff.getlijsten()
    gffb = gffb_bestand(bestandsnaam_gbff)
    gffb_gene, gffb_cds, gffb_rrna = gffb.getdiconaries()

    #sequentie_genoom = sequentie_uit_bestand(sequentie_bestand)
    #lijst_met_gene_cds = informatie_uit_gene_cds_bestand(gene_cds_bestand)
    #gene, cds = maken_dictonary(lijst_met_gene_cds)

    #grafieken(inhoud_gff,sequentie_genoom,cds)
    #gevonden_patronen = patronen(cds)
    p = Patronen(gff_cds,gffb_gene,gffb_cds)
    patronen = p.getpatronen()

    g = gui_maken(patronen,gff_totale_inhoud,gffb_cds)

main()
