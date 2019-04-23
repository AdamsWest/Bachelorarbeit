Leistungsuntersuchung von elektrisch, propellergetriebenen Flugger�ten

Alle Untersuchungen von Leistungsparametern sind �hnlich aufgebaut. Bei einigen Untersuchungen sind zus�tzliche Skripte vorhanden, um gesonderte Zusammenh�nge zu untersuchen.


- die Leistungsberechnung wird mit dem Startskript-Leistungsberechnung gestartet
- Startskript_Leistungsberechnung: Definition der Flugger�te, der Motor, der Propeller, der Batterie und der Umgebungsparameter
- Leistungsberechnung: Hauptprogramm, in dem die Leistungsberechnung durchgef�hrt wird, nicht ver�ndern!
- Batterie: Funktion, in der der Batteriezustand berechnet wird, nicht ver�ndern!
- Batterie_parameter: Funktion, Berechnung wichtiger Batterieparameter zur Bestimmung der Batteriespannung, nicht ver�ndern!
- ESC: Funktion, Berechnung des Motorreglerzustands, nicht ver�ndern!
- FlaechenflugzeugAerodynamik: Funktion, ermittelt den erforderlichen Schub f�r ein Fl�chenflugzeug, nicht ver�ndern!
- Motor: Funktion, Berechnung des Motorzustands, nicht ver�ndern!
- Motordata: Funktion, Entnahme aller Motorkennparameter aus der AXI-Datenbank, nicht ver�ndern!
- MulticopterAerodynamik: Funktion, ermittelt den erforderlichen Schub f�r einen Multicopter, nicht ver�ndern!
- Normcell: Funktion, Berechnung einer normierten Batteriezelle aus Batteriedatenbank, nicht ver�ndern!
- Propeller: Funktion, Interpolation der Propellerdrehzahl und des Propellerdrehmoments aus den Propellerkennfeldern, nicht ver�ndern!
- Propeller_map: Funktion, Transformation der APC-Propellerkennfelder in �quidistante Geschwindigkeitsabst�nde, nicht ver�ndern!


Datenbanken:
- axi_motor_db: Motordatenbank
- DATA_APC: Propellerkennfelddatenbank
- Elektromodellflug: Batteriedatenbank


Eine ausf�hrliche Dokumentation der Funktionen und des Programmablaufs wird in ForceCalc vorgenommen.