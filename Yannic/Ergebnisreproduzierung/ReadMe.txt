Leistungsuntersuchung von elektrisch, propellergetriebenen Fluggeräten

Alle Untersuchungen von Leistungsparametern sind ähnlich aufgebaut. Bei einigen Untersuchungen sind zusätzliche Skripte vorhanden, um gesonderte Zusammenhänge zu untersuchen.


- die Leistungsberechnung wird mit dem Startskript-Leistungsberechnung gestartet
- Startskript_Leistungsberechnung: Definition der Fluggeräte, der Motor, der Propeller, der Batterie und der Umgebungsparameter
- Leistungsberechnung: Hauptprogramm, in dem die Leistungsberechnung durchgeführt wird, nicht verändern!
- Batterie: Funktion, in der der Batteriezustand berechnet wird, nicht verändern!
- Batterie_parameter: Funktion, Berechnung wichtiger Batterieparameter zur Bestimmung der Batteriespannung, nicht verändern!
- ESC: Funktion, Berechnung des Motorreglerzustands, nicht verändern!
- FlaechenflugzeugAerodynamik: Funktion, ermittelt den erforderlichen Schub für ein Flächenflugzeug, nicht verändern!
- Motor: Funktion, Berechnung des Motorzustands, nicht verändern!
- Motordata: Funktion, Entnahme aller Motorkennparameter aus der AXI-Datenbank, nicht verändern!
- MulticopterAerodynamik: Funktion, ermittelt den erforderlichen Schub für einen Multicopter, nicht verändern!
- Normcell: Funktion, Berechnung einer normierten Batteriezelle aus Batteriedatenbank, nicht verändern!
- Propeller: Funktion, Interpolation der Propellerdrehzahl und des Propellerdrehmoments aus den Propellerkennfeldern, nicht verändern!
- Propeller_map: Funktion, Transformation der APC-Propellerkennfelder in äquidistante Geschwindigkeitsabstände, nicht verändern!


Datenbanken:
- axi_motor_db: Motordatenbank
- DATA_APC: Propellerkennfelddatenbank
- Elektromodellflug: Batteriedatenbank


Eine ausführliche Dokumentation der Funktionen und des Programmablaufs wird in ForceCalc vorgenommen.