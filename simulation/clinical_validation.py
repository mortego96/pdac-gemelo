import sys
import os
import numpy as np
import time

# Permite ejecutar desde cualquier lado
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from simulation.tumor_model import TumorModel


class ClinicalValidationSuite:
    def __init__(self, cohort_size=3, sim_hours=1440):
        # 1440 horas = 60 días in silico (suficiente para evaluar TGI inicial y resistencias)
        self.cohort_size = cohort_size
        self.sim_hours = sim_hours
        self.results = {}
        
    def _run_virtual_patient(self, protocol_name, mut_profile=None):
        """Corre la simulación para 1 paciente virtual y extrae el TGI (Tumor Growth Inhibition)."""
        m = TumorModel(width=50, height=50, n_cancer=300, n_caf=150)
        
        # Forzar mutaciones necesarias para el caso de uso
        if mut_profile:
            for cell in m.agents:
                if hasattr(cell, 'cell_type') and cell.cell_type == 'cancer':
                    cell.signaling.mutations.update(mut_profile)
                    
        m.set_protocol(protocol_name)
        
        initial_tumor_burden = len([a for a in m.agents if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive])
        
        # Ejecutar simulación sin output visual
        for _ in range(self.sim_hours):
            m.step()
            
        final_tumor_burden = len([a for a in m.agents if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive])
        
        # Calcular RECIST proxy (variación %)
        tgi_pct = ((final_tumor_burden - initial_tumor_burden) / initial_tumor_burden) * 100
        
        # Extraer biomarcadores promedio para revisar resistencias
        cancer_cells = [a for a in m.agents if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive]
        biomarkers = {}
        if cancer_cells:
            # Aggregate nodes
            for key in cancer_cells[0].signaling.nodes.keys():
                vals = [c.signaling.nodes[key] for c in cancer_cells]
                biomarkers[key] = np.mean(vals)
                
        return tgi_pct, biomarkers

    def run_arm(self, arm_name, protocol, expected_tgi_range, mut_profile=None):
        """Ejecuta un brazo entero de la prueba clínica (N pacientes)."""
        print(f"\\n➤ Iniciando Brazo Clínico: {arm_name} (Protocolo: {protocol})")
        print(f"  Pacientes Virtuales (N={self.cohort_size}), Duración={self.sim_hours//24} días. Esperado: {expected_tgi_range[0]}% a {expected_tgi_range[1]}%")
        
        arm_tgis = []
        arm_biomarkers = []
        
        start_t = time.time()
        for i in range(self.cohort_size):
            tgi, bios = self._run_virtual_patient(protocol, mut_profile)
            arm_tgis.append(tgi)
            arm_biomarkers.append(bios)
            print(f"    Paciente {i+1}: Variación tumoral = {tgi:+.1f}%")
            
        avg_tgi = np.mean(arm_tgis)
        std_tgi = np.std(arm_tgis)
        
        passed = expected_tgi_range[0] <= avg_tgi <= expected_tgi_range[1]
        
        print(f"  [RESULTADO] Promedio TGI: {avg_tgi:+.1f}% ± {std_tgi:.1f}%")
        if passed:
            print("  ✅ VERIFICADO: El gemelo digital reproduce fielmente la cinética observada en ensayos clínicos reales.")
        else:
            print("  ❌ FALLO DE CALIBRACIÓN: La respuesta in silico está fuera del rango empírico (falso positivo/negativo).")
            
        self.results[arm_name] = {
            'avg_tgi': avg_tgi,
            'passed': passed,
            'biomarkers': arm_biomarkers[0] if arm_biomarkers else {}
        }

    def validate_resistances(self):
        print("\\n➤ Verificando Precisión de Vías de Resistencia Moleculares (In Silico vs 2024 Literature)")
        
        # Revisar FOLFIRINOX -> NOTCH/DDR
        if 'FOLFIRINOX' in self.results:
            bios = self.results['FOLFIRINOX']['biomarkers']
            if bios.get('NOTCH_active', 0) > 0.4 and bios.get('DDR_active', 0) > 0.4:
                print("  ✅ PASS | FOLFIRINOX indujo resistencia por vía NOTCH/DDR (GALNT5-axis).")
            else:
                print(f"  ❌ FAIL | FOLFIRINOX falló al activar NOTCH/DDR (NOTCH={bios.get('NOTCH_active',0):.2f}, DDR={bios.get('DDR_active',0):.2f})")
                
        # Revisar KRASi -> Macropinocitosis y MYC
        if 'KRASi (RMC-6236)' in self.results:
            bios = self.results['KRASi (RMC-6236)']['biomarkers']
            if bios.get('macropinocytosis', 0) > 0.45 and bios.get('MYC_active', 0) > 0.35:
                print("  ✅ PASS | KRASi indujo escape metabólico vía Macropinocitosis y amplificación MYC (AGER-DIAPH1 axis).")
            else:
                print(f"  ❌ FAIL | KRASi falló escape metabólico (Macro={bios.get('macropinocytosis',0):.2f}, MYC={bios.get('MYC_active',0):.2f})")

        # Revisar Olaparib -> OXPHOS y BRCA
        if 'Olaparib' in self.results:
             bios = self.results['Olaparib']['biomarkers']
             if bios.get('OXPHOS_active', 0) > 0.35 and bios.get('BRCA_functional', 0) > 0.4:
                 print("  ✅ PASS | Olaparib indujo shift a OXPHOS y reversión fenotípica de BRCA (Restauración HRD).")
             else:
                 print(f"  ❌ FAIL | PARPi falló recableado OXPHOS (OXPHOS={bios.get('OXPHOS_active',0):.2f}, BRCA_func={bios.get('BRCA_functional',0):.2f})")


if __name__ == '__main__':
    print("==================================================================================")
    print(" SUITE DE VALIDACIÓN CLÍNICA EXTERNA — GEMELO DIGITAL PDAC (vs Ensayos 2024-2025)")
    print("==================================================================================\\n")
    
    suite = ClinicalValidationSuite(cohort_size=3, sim_hours=1440) # 2 meses in silico aprox
    
    # BRAZO 1: Historia Natural (Control no tratado) -> Esperamos crecimiento masivo (>100%)
    suite.run_arm("Control (Historia Natural)", "Sin_tratamiento", expected_tgi_range=(80, 500))
    
    # BRAZO 2: FOLFIRINOX (PRODIGE 4) -> Esperamos ORR 31%. Reducción o estabilidad duradera (-60% a +20%).
    suite.run_arm("FOLFIRINOX", "FOLFIRINOX", expected_tgi_range=(-60, 20))
    
    # BRAZO 3: Gem/NabPac (MPACT) -> Esperamos ORR 23%. Reducción moderada o enfermedad estable.
    suite.run_arm("Gemcitabina + Nab-Paclitaxel", "GemNabPac", expected_tgi_range=(-45, 30))
    
    # BRAZO 4: Inhibidor KRAS G12D/Pan-KRAS -> Fuerte "tumor shrinkage" inicial pero con escape.
    suite.run_arm("KRASi (RMC-6236)", "KRASi_mono", expected_tgi_range=(-70, 10), mut_profile={'KRAS': 'HET_MUT', 'KRAS_variant': 'G12V'})

    # BRAZO 5: PARP Inhibitor en BRCA mutado (POLO) -> Mantenimiento, no suele reducir agresivamente pero bloquea crecimiento.
    suite.run_arm("Olaparib", "Olaparib_maint", expected_tgi_range=(-30, 10), mut_profile={'BRCA': 'HOM_LOSS'})

    suite.validate_resistances()
    print("\\n==================================================================================")
    print(" FINALIZADA LA AUDITORÍA DE VALIDACIÓN. Revisa los resultados para calibrar si es necesario.")
    print("==================================================================================")
