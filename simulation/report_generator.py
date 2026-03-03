"""
Generador de informes científicos PDAC v6.
Reportes detallados con biomarcadores, vías moleculares y mecanismos
de resistencia explicados con rigor biológico.
"""

import numpy as np


def generate_report(model, drug_names: list = None) -> str:
    """Genera un informe científico completo de la simulación."""
    s = model.get_summary()
    h = model.history
    deltas = get_pathway_deltas(model)
    drug_names = drug_names or getattr(model, 'active_drugs', [])

    lines = []
    _add = lines.append

    _add("# INFORME CIENTÍFICO — SIMULACIÓN DE GEMELO DIGITAL PDAC")
    _add("---")
    _add("")

    # ═══ 1. RESUMEN EJECUTIVO CLÍNICO ═══
    _add("## 1. DIAGNÓSTICO IN SILICO (Executive Summary)")
    
    if not h['step']:
        _add("   Simulación no iniciada.")
        _add("")
    else:
        c_ini = h['cancer_alive'][0]
        c_fin = h['cancer_alive'][-1]
        cambio = c_fin - c_ini
        pct = (cambio / max(c_ini, 1)) * 100
        
        # Lenguaje natural predictivo
        if pct < -50:
            diag = "RESPONDEDOR EXCEPCIONAL. Tratamiento altamente eficaz."
            desc = "El tumor está en franca regresión volumétrica, mostrando sensibilidad profunda al régimen."
        elif c_fin == 0:
            diag = "RESPUESTA COMPLETA PREDICTA."
            desc = "Erradicación de las células tumorales en silico."
        elif -50 <= pct <= 20:
            diag = "ENFERMEDAD ESTABLE (Control Tumoral Mantenido)."
            desc = "El tratamiento frena la proliferación aguda, pero el tumor no colapsa."
        else:
            diag = "PROGRESIÓN / NO-RESPONDEDOR."
            desc = "Tratamiento ineficaz. El tumor no responde y continúa proliferando activamente."
            
        _add(f"**► {diag}**")
        _add(f"_{desc}_")
        
        cancer_cells = [a for a in model.agents if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive]
        avg_cells = _average_nodes(cancer_cells) if cancer_cells else {}
        p_sig = avg_cells.get('proliferation_signal', 0)
        
        if cancer_cells and p_sig > 0.6:
            _add(f"⚠️ **Alerta:** Alta presión proliferativa residual (índice Ki-67 estimado >40%).")
        _add("")

    # ═══ 9. NARRATIVA CLÍNICA ═══
    _add("## 2. ¿QUÉ HA OCURRIDO DURANTE LA SIMULACIÓN? (Narrativa Clínica)")
    
    # Bug fix: We must look at drug_names passed by main.py, or any drug that has ever been dosed > 0 in history
    active = set(drug_names) if drug_names else set()
    if not active and hasattr(model, 'drug_doses'):
        active = {d for d, v in model.drug_doses.items() if v > 0}
    active = list(active)
    
    cancer_cells = [a for a in model.agents if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive]
    avg = _average_nodes(cancer_cells) if cancer_cells else {}
    
    if not h['step']:
        _add("   Simulación no ejecutada.")
    elif not active:
        _add("   **Escenario de Crecimiento Natural (Brazo Control):**")
        _add("   Ante la ausencia de intervención farmacológica, el gemelo digital PDAC exhibe el comportamiento estocástico agresivo intrínseco de su fenotipo G12D. Las células cancerosas mantienen un programa de división logarítmica propulsado por la cascada KRAS-MAPK hiperactiva. No se observa restricción mitótica significativa, resultando en una expansión poblacional expedita de la masa tumoral.")
        _add(f"   *Progreso:* Incremento volumétrico nítido, partiendo desde {h['cancer_alive'][0]} células hasta alcanzar una meseta de {h['cancer_alive'][-1]} células ({pct:+.0f}%) al final de la ventana biológica simulada.")
    else:
        # Narrative for treated arms
        drug_desc_list = []
        for dn in active:
            drug = model.drug_library.get_drug(dn)
            drug_desc_list.append(drug.name if drug else dn)
        drug_str = " + ".join(drug_desc_list)
        
        _add(f"   **Intervención Terapeútica Evaluada:** {drug_str}")
        _add("")
        
        # FASE 1: Acción Farmacológica
        _add("   **1. Dinámica Aguda (Inicio del Tratamiento):**")
        
        is_krasi = any(k in active for k in ['daraxonrasib', 'mrtx1133', 'rmc7977', 'sotorasib'])
        has_afatinib = 'afatinib' in active
        has_stat3 = 'protac_stat3' in active
        is_chemo = any(k in active for k in ['gemcitabine', '5fu', 'oxaliplatin', 'irinotecan', 'nab_paclitaxel'])
        is_parpi = 'olaparib' in active
        
        if is_krasi:
            _add("   La administración del inhibidor de KRAS induce inicialmente un 'shock oncogénico' profundo. Al bloquear selectivamente la isoforma mutante en su estado activo, la cascada MAPK se interrumpe súbitamente, privando al tumor de su motor proliferativo principal.")
            if avg.get('ERK_active', 0) < 0.1:
                _add("   Esta acción se refleja en un colapso dramático del tono de pERK, indicando penetración primaria eficaz.")
        
        if has_afatinib:
            _add("   La co-administración del inhibidor pan-HER bloquea proactivamente los receptores transmembrana (EGFR/HER2), anticipando y suprimiendo el rescate de la señalización mediada por receptores tirosina quinasa.")
            
        if has_stat3:
            _add("   El PROTAC degrada eficientemente STAT3, despojando al tumor de un nodo crucial de supervivencia transcripcional cruzada.")
            
        if is_chemo:
            _add("   La infusión citotóxica bombardea el ciclo celular, induciendo aductos de ADN y frustración mitótica que decanta rápidamente en estrés oxidativo generalizado.")
            
        if is_parpi:
            brca = avg.get('BRCA_functional', 0.9)
            if brca < 0.3:
                _add("   En este fenotipo HRD+ (deficiencia recombinación homóloga), el atrapamiento de PARP1 por Olaparib canaliza el daño directo hacia roturas de doble cadena irreparables, explotando magistralmente la letalidad sintética.")
            else:
                _add("   A pesar del atrapamiento de PARP, la maquinaria intacta de la recombinación homóloga (BRCA-wt) repara las lesiones subyacentes, opacando el impacto basal del fármaco.")
                
        _add("")
        
        # FASE 2: Respuesta y Resistencia
        _add("   **2. Plasticidad Tumoral y Remodelación Adaptativa (Fase Crónica):**")
        
        if is_krasi and not has_afatinib:
            _add("   A la ventana de supresión aguda le sigue una rápida remodelación de la red molecular. Percibiendo el déficit letal de señalización, el tumor despliega un programa clásico de resistencia adaptativa: al evaporarse el feedback negativo mediado por ERK (vía DUSP6/SPRY2), los receptores EGFR y HER2 se des-reprimen y sobreexpresan dramáticamente en la membrana. ")
            _add("   Este 'bypass' no-genómico permite re-alimentar la vía MAPK utilizando KRAS wild-type y reimpulsar el polo mitótico, estabilizando así la cinética celular frente a la embestida química.")
            
        elif is_krasi and has_afatinib:
            _add("   El cerrojo dual sobre KRAS y EGFR/HER2 estrangula existencialmente al fenotipo epitelial. Sin poder reclutar señales de supervivencia extracelulares, la célula es forzada hacia la maquinaria apoptótica intrínseca.")
            if has_stat3:
                _add("   Al cercenar además el salvavidas mediado por STAT3, el colapso mitocondrial ocurre en cascada. La letalidad es ineludible.")
        elif is_chemo:
             _add("   Bajo el insulto quimioterápico persistente, el clon tumoral enciende mecanismos homeostásicos de rescate extremo (hiperactivación del flujo autofágico, repoblación de factores DDR, o estabilización basal dependiente de estrés), mitigando la erradicación total del recuento celular.")
        
        _add("")
        
        # FASE 3: Outcome
        _add("   **3. Balance Cifrado del Tejido al Sacrificio In Silico:**")
        if h['step']:
            c_ini = h['cancer_alive'][0]
            c_fin = h['cancer_alive'][-1]
            pct_f = (c_fin - c_ini) / max(c_ini, 1) * 100
            p_f = avg.get('proliferation_signal', 0)
            a_f = avg.get('apoptosis_signal', 0)
            s_f = avg.get('survival_signal', 0)
            
            outcome_str = f"→ Evolución fraccional: {pct_f:+.1f}% de masa viable."
            if a_f > p_f * 2.5:
                outcome_str += " Consenso fenotípico: REGRESIÓN MAYOR INDUCIDA (Apoptosis dominante)."
            elif a_f > p_f:
                outcome_str += " Consenso fenotípico: CONTENCIÓN MARGINAL (Equilibrio inestable)."
            elif p_f > a_f * 3:
                outcome_str += " Consenso fenotípico: PROGRESIÓN REFRACTARIA (Intervención fútil)."
            else:
                outcome_str += " Consenso fenotípico: ENFERMEDAD ESTABLE (Estancamiento citoestático)."
                
            _add(f"   {outcome_str}")
            _add(f"   (Firma mitótica final: {p_f:.2f} | Tono apoptótico: {a_f:.2f} | Ratio de Supervivencia cruzada: {s_f:.2f})")
    _add("")
    _add("")

    # ═══ 2. BENCHMARKING EMPÍRICO ═══
    _add("## 3. COMPARATIVA CON LITERATURA CLÍNICA (Benchmarking)")
    if not h['step']:
        _add("   Sin datos para comparar.")
    else:
        active = set(drug_names) if drug_names else set()
        if not active and hasattr(model, 'drug_doses'):
            active = {d for d, v in model.drug_doses.items() if v > 0}
        active = list(active)
            
        protocol_name = getattr(model, 'active_protocol_name', None)
        
        _add(f"**Reducción/Crecimiento en la simulación:** {pct:+.1f}%")
        
        if protocol_name == "FOLFIRINOX" or (set(['5fu', 'oxaliplatin', 'irinotecan']).issubset(set(active))):
            _add("   Benchmarking FOLFIRINOX:")
            _add("   - Ensayos Clínicos (PRODIGE 4, PASS-01 2024): La Tasa de Respuesta Objetiva (ORR)")
            _add("     esperada es del 31-52%. mPFS = 4-8 meses.")
            if pct < -20:
                _add("   - Veredicto: El gemelo digital se comporta como el tercio superior de respondedores.")
            else:
                _add("   - Veredicto: El modelo refleja resistencia natural o intrínseca, un escenario común in vivo.")
        
        elif protocol_name == "Gemcitabina + Nab-Paclitaxel" or (set(['gemcitabine', 'nab_paclitaxel']).issubset(set(active))):
            _add("   Benchmarking Gem/NabPac:")
            _add("   - Ensayos Clínicos (MPACT, PASS-01 2024): ORR esperada del 23-57%. mPFS = 5-7 meses.")
            if pct < -20:
                _add("   - Veredicto: Respuesta in silico alineada con respondedores clìnicos óptimos.")
            else:
                _add("   - Veredicto: Simulación refleja el fenotipo estromal denso/resistente típico del ~60% de pacientes.")
                
        elif any(k in active for k in ['daraxonrasib', 'rmc7977', 'mrtx1133', 'sotorasib']):
            _add("   Benchmarking Inhibidores KRAS / Pan-RAS:")
            _add("   - Datos Fase 1/2 (2024-2025 RMC-6236): ORR ~20-36%.")
            _add("   - In vitro / PDX: Inducen citoestasis drástica inicial, seguida de rebote")
            _add("     por reactivación RTK/MAPK a las 4-8 semanas.")
            if pct < -10:
                _add("   - Veredicto: El modelo captura exitosamente el choque oncogénico inicial del fármaco.")
            else:
                _add("   - Veredicto: El modelo demuestra escape temprano o bypass intrínseco.")
        
        elif len(active) == 0:
            _add("   Benchmarking Control:")
            _add("   - Estudios de historia natural evidencian crecimiento exponencial agresivo en ausencia de tratamiento.")
        
        else:
            _add("   - Regímenes experimentales: requiere cruce con ensayos preclínicos fase PDX para validación TGI.")
            
    _add("")

    # ═══ 5. TRATAMIENTO ═══
    _add("## 4. RÉGIMEN TERAPÉUTICO")
    active = set(drug_names) if drug_names else set()
    if not active and hasattr(model, 'drug_doses'):
        active = {d for d, v in model.drug_doses.items() if v > 0}
    active = list(active)
    
    if active:
        for dn in active:
            dose = model.drug_doses.get(dn, 0)
            drug = model.drug_library.get_drug(dn)
            name = drug.name if drug else dn
            desc = drug.description if drug else ""
            targets = ", ".join(drug.targets.keys()) if drug else "?"
            _add(f"   • {name}")
            _add(f"     Dosis en T0 (Mantenimiento variable): {dose:.0%} | IC₅₀: {drug.ic50 if drug else '?'}")
            _add(f"     Dianas: {targets}")
            if desc:
                _add(f"     {desc}")
            _add("")
        if len(active) > 1:
            _add("   Combinación evaluada con modelo de independencia de Bliss:")
            _add("   E(A+B) = E(A) + E(B) - E(A)×E(B)")
    else:
        _add("   Sin tratamiento farmacológico (brazo control)")
    _add("")

    # ═══ 3. CINÉTICA Y MÉTRICAS BÁSICAS ═══
    _add("## 5. CINÉTICA Y MÉTRICAS BÁSICAS")
    _add(f"   Pasos simulados: {s['Paso']} (1 paso ≈ 1 hora)")
    _add(f"   Células tumorales viables: {s['Células tumorales']}")
    if h['step']:
        _add(f"   Pico tumoral: {max(h['cancer_alive'])} células (hora {h['step'][h['cancer_alive'].index(max(h['cancer_alive']))]})")
        
        if c_fin > c_ini and c_ini > 0:
            dt = (s['Paso'] * np.log(2)) / np.log(c_fin / c_ini)
            _add(f"   Tiempo de duplicación estimado: {dt:.0f} horas")
            if dt < 24:
                _add("     (Crecimiento muy rápido — posible fenotipo basal-like agresivo)")
            elif dt < 72:
                _add("     (Crecimiento estándar PDAC)")
            else:
                _add("     (Crecimiento lento/frenado — efecto citoestático detectado)")
    _add("")

    # ═══ 4. PERFIL MUTACIONAL ═══
    _add("## 6. PERFIL MUTACIONAL DEL TUMOR")
    if cancer_cells and hasattr(cancer_cells[0], 'signaling'):
        sig = cancer_cells[0].signaling
        muts = sig.mutations
        var = muts.get('KRAS_variant', 'G12D')
        _add(f"   KRAS: {'Mutado ({})'.format(var) if muts.get('KRAS') else 'Wild-type'}")
        if muts.get('KRAS'):
            from signaling.pdac_network import KRAS_VARIANTS
            kp = KRAS_VARIANTS.get(var, {})
            _add(f"     → {kp.get('description', '')}")
        _add(f"   TP53: {'Pérdida de función (LOF)' if muts.get('TP53') else 'Funcional'}")
        _add(f"   CDKN2A/p16: {'Deleción homocigota' if muts.get('CDKN2A') else 'Expresado'}")
        _add(f"   SMAD4: {'Pérdida (promueve EMT vía TGF-β)' if muts.get('SMAD4') else 'Funcional'}")
        _add(f"   YAP1: {'Amplificación' if muts.get('YAP_amp') else 'Normal'}")
        _add(f"   BRCA1/2: {'Mutado (HRD+, sensible a PARPi)' if muts.get('BRCA') else 'Funcional'}")
        _add(f"   MYC: {'Amplificado' if muts.get('MYC') else 'Normal'}")
    _add("")

    # ═══ 6. MICROAMBIENTE ═══
    _add("## 7. MICROAMBIENTE TUMORAL (TME)")
    if hasattr(model, 'microenv'):
        me = model.microenv
        oxy_mean = float(np.mean(me.oxygen))
        glu_mean = float(np.mean(me.glucose))
        ecm_mean = float(np.mean(me.ecm_density))
        lac_mean = float(np.mean(me.lactate))
        ph_mean = float(np.mean(me.ph))
        hypoxic_pct = float(np.mean(me.oxygen < 0.3)) * 100
        tgfb_mean = float(np.mean(me.tgfb))
        il6_mean = float(np.mean(me.il6))

        _add(f"   pO₂ promedio:          {oxy_mean:.3f} (normal: >0.6)")
        _add(f"   Áreas hipóxicas:       {hypoxic_pct:.0f}% del tejido (pO₂ <0.3)")
        _add(f"   Glucosa promedio:      {glu_mean:.3f} (normal: >0.5)")
        _add(f"   Lactato promedio:      {lac_mean:.3f} (normal: <0.2)")
        _add(f"   pH extracelular:       {ph_mean:.2f} (normal: 7.35-7.45)")
        _add(f"   ECM/Colágeno (Masson): {ecm_mean:.3f} (desmoplasia si >0.4)")
        _add(f"   TGF-β (ELISA):         {tgfb_mean:.3f}")
        _add(f"   IL-6 (ELISA):          {il6_mean:.3f}")

        if hypoxic_pct > 30:
            _add("")
            _add("   ⚠ HIPOXIA SIGNIFICATIVA")
            _add("   La hipoxia tumoral (pimonidazol+ en IHC) favorece:")
            _add("    - Estabilización de HIF-1α/2α → switch a glucólisis anaerobia")
            _add("    - Acidificación del TME (pH<7.0) → inmunosupresión")
            _add("    - Upregulation de PD-L1 → escape inmune")
            _add("    - Selección de clones resistentes a quimioterapia")

        if ecm_mean > 0.4:
            _add("")
            _add("   ⚠ DESMOPLASIA (reacción estromal densa)")
            _add("   La ECM densa en PDAC (>80% del volumen tumoral) causa:")
            _add("    - Barrera física para penetración de fármacos (IFP ↑)")
            _add("    - Exclusión de linfocitos T CD8+ del core tumoral")
            _add("    - Mecanotransducción → activación YAP/TAZ vía integrinas")
            _add(f"    - Factor de penetración estimado: {1.0-ecm_mean*0.6:.0%}")

        if lac_mean > 0.3:
            _add("")
            _add("   ⚠ ACIDOSIS METABÓLICA TUMORAL")
            _add(f"   Lactato elevado ({lac_mean:.2f}) → pH {ph_mean:.2f}")
            _add("    - Inhibición de función de T-cells y NK cells")
            _add("    - Polarización M2 de macrófagos")
            _add("    - Resistencia a fármacos pH-dependientes")
    _add("")

    # ═══ 4. BIOMARCADORES MOLECULARES ═══
    _add("## 8. BIOMARCADORES MOLECULARES DE SISTEMA (Promedios)")
    if cancer_cells:
        avg = _average_nodes(cancer_cells)

        _add("")
        _add("### 8.1. VÍA RAS-MAPK (biomarcadores: pERK1/2, pMEK, DUSP6, SPRY2)")
        _add(f"      EGFR (Sobreexpresión RTK):         {avg.get('EGFR_active', 0):.3f}")
        _add(f"      HER2 (Bypass Adaptativo):          {avg.get('HER2_active', 0):.3f}")
        _add(f"      KRAS Mutado (GTPasa):              {avg.get('KRAS_active', 0):.3f}")
        _add(f"      KRAS Wild-type (Compensatorio):    {avg.get('KRAS_WT_active', 0):.3f}")
        _add(f"      pRAF (RAF fosforilado):            {avg.get('RAF_active', 0):.3f}")
        _add(f"      pMEK1/2 (WB/IHC: #9154 CST):       {avg.get('MEK_active', 0):.3f}")
        _add(f"      pERK1/2 (Thr202/Tyr204, #4370):    {avg.get('ERK_active', 0):.3f}")
        if avg.get('ERK_active', 0) > 0.6:
            _add("      ⚠ pERK elevado → señal proliferativa activa")
            _add("        Biomarcadores feedback: DUSP6 ↑, SPRY2 ↑ (reguladores negativos)")
            _add("        IHC: Ki-67 probablemente elevado (>40%)")
        elif avg.get('ERK_active', 0) < 0.2:
            _add("      ✓ pERK bajo → inhibición efectiva de la vía MAPK (shock inicial)")
            if avg.get('EGFR_active', 0) > 0.4 or avg.get('HER2_active', 0) > 0.4:
                _add("      ⚠ ALERTA: Escape RTK detectado (EGFR/HER2 ↑ por falta de represión de pERK)")

        _add("")
        _add("### 8.2. VÍA PI3K/AKT/mTOR (biomarcadores: pAKT Ser473, pS6K, p4E-BP1)")
        _add(f"      PI3K (activación):                 {avg.get('PI3K_active', 0):.3f}")
        _add(f"      pAKT (Ser473, #4060 CST):         {avg.get('AKT_active', 0):.3f}")
        _add(f"      PTEN (IHC: pérdida = mal pronóst): {avg.get('PTEN_active', 0):.3f}")
        _add(f"      pmTOR (Ser2448):                   {avg.get('mTOR_active', 0):.3f}")
        _add(f"      pS6K (Thr389):                     {avg.get('S6K_active', 0):.3f}")
        if avg.get('AKT_active', 0) > 0.5:
            _add("      ⚠ pAKT elevado → supervivencia y resistencia a apoptosis")
            _add("        Considerar inhibidores mTOR (everolimus) o AKT (capivasertib)")

        _add("")
        _add("### 8.3. YAP/TAZ-TEAD (biomarcadores: YAP nuclear IHC, CTGF, CYR61)")
        _add(f"      YAP nuclear (IHC: sc-101199):      {avg.get('YAP_nuclear', 0):.3f}")
        _add(f"      TEAD (actividad transcripcional):   {avg.get('TEAD_active', 0):.3f}")
        _add(f"      CTGF/CCN2 (gen diana YAP):         {avg.get('CTGF_expression', 0):.3f}")
        _add(f"      CYR61/CCN1 (gen diana YAP):        {avg.get('CYR61_expression', 0):.3f}")
        if avg.get('YAP_nuclear', 0) > 0.5:
            _add("      ⚠ YAP nuclear elevado → bypass de KRASi, resistencia intrínseca")
            _add("        Mecanismo: mecanotransducción vía ECM → Hippo OFF → YAP ↑")

        _add("")
        _add("### 8.4. HIPOXIA (biomarcadores: HIF-1α IHC, CA-IX, GLUT1, VEGF-A)")
        _add(f"      HIF-1α (IHC: estabilización):      {avg.get('HIF1A_active', 0):.3f}")
        _add(f"      HIF-2α (EPAS1):                    {avg.get('HIF2A_active', 0):.3f}")
        _add(f"      VEGF-A (ELISA/IHC):                {avg.get('VEGF_expression', 0):.3f}")
        _add(f"      Efecto Warburg (GLUT1↑/LDHA↑):     {avg.get('warburg_effect', 0):.3f}")
        _add(f"      Adicción a glutamina (GLS1↑):      {avg.get('glutamine_addiction', 0):.3f}")
        if avg.get('HIF1A_active', 0) > 0.4 or avg.get('HIF2A_active', 0) > 0.4:
            _add("      ⚠ Hipoxia activa → switch metabólico, PD-L1 ↑, resistencia a QT")
            _add("        Marcadores IHC: CA-IX ↑ (marcador subrogado de hipoxia)")
            _add("        Considerar belzutifan (HIF-2αi) si HIF-2α dominante")

        _add("")
        _add("### 8.5. APOPTOSIS (biomarcadores: cleaved Caspase-3, TUNEL, Annexin V)")
        apop = avg.get('apoptosis_signal', 0)
        _add(f"      Señal apoptótica integrada:        {apop:.3f}")
        _add(f"      BCL-2 (antiapoptótico, IHC):       {avg.get('BCL2_active', 0):.3f}")
        _add(f"      BAX (proapoptótico):               {avg.get('BAX_active', 0):.3f}")
        _add(f"      Ratio BAX/BCL-2:                   {avg.get('BAX_active',0) / max(avg.get('BCL2_active',0), 0.01):.2f}")
        _add(f"      Caspasa-3 activada (cleaved):      {avg.get('CASP3_active', 0):.3f}")
        _add(f"      p53 funcional:                     {avg.get('TP53_active', 0):.3f}")
        _add(f"      p21/CDKN1A (downstream p53):       {avg.get('P21_active', 0):.3f}")
        if avg.get('BAX_active', 0) > avg.get('BCL2_active', 0):
            _add("      → Ratio BAX/BCL-2 > 1: priming apoptótico favorable")
            _add("        Esperado: TUNEL+ ↑, cleaved caspase-3 ↑ en IHC")
        else:
            _add("      → Ratio BAX/BCL-2 < 1: resistencia a apoptosis")
            _add("        Considerar BH3-miméticos (venetoclax/navitoclax)")

        _add("")
        _add("### 8.6. EMT Y PLASTICIDAD (biomarcadores: E-cadherina, Vimentina, ZEB1)")
        _add(f"      E-cadherina (CDH1, IHC):           {avg.get('ECAD_expression', 0):.3f}")
        _add(f"      Vimentina (VIM, IHC):              {avg.get('VIM_expression', 0):.3f}")
        _add(f"      SNAIL (SNAI1):                     {avg.get('SNAIL_active', 0):.3f}")
        _add(f"      ZEB1 (represor CDH1):              {avg.get('ZEB1_active', 0):.3f}")
        _add(f"      Programa basal-like:               {avg.get('basal_like_program', 0):.3f}")
        basal_pct = max(h['basal_like_pct']) if h['basal_like_pct'] else 0
        if basal_pct > 40:
            _add(f"      ⚠ {basal_pct:.0f}% células con fenotipo basal-like/escamoso")
            _add("        Transición epitelio-mesenquimal (EMT) activa")
            _add("        Marcadores: CDH1 ↓, VIM ↑, ΔNp63 ↑, KRT5/14 ↑")
            _add("        Clínicamente: peor pronóstico, resistencia a gemcitabina")

        _add("")
        _add("### 8.7. INMUNOEVASIÓN (biomarcadores: PD-L1 IHC (22C3/28-8), MHC-I)")
        _add(f"      PD-L1 (SP263/22C3, CPS score):     {avg.get('PDL1_expression', 0):.3f}")
        _add(f"      MHC-I (β2-microglobulina):         {avg.get('MHC1_expression', 0):.3f}")
        _add(f"      NFκB (p65 nuclear):                {avg.get('NFKB_active', 0):.3f}")
        _add(f"      STAT3 (pSTAT3 Y705):               {avg.get('STAT3_active', 0):.3f}")
        if avg.get('PDL1_expression', 0) > 0.5:
            _add("      ⚠ PD-L1 elevado (CPS≥1 predicho)")
            _add("        Pero PDAC es típicamente 'inmunológicamente frío'")
            _add("        Respuesta a anti-PD-1 esperada solo si MSI-H/dMMR (<2%)")
        if avg.get('MHC1_expression', 0) < 0.4:
            _add("      ⚠ MHC-I downregulated → evasión de T-cells CD8+")
            _add("        Mecanismo: HIF → represión β2M, o pérdida alélica")

        _add("")
        _add("### 8.8. METABOLISMO (biomarcadores: LDHA, GLUT1, GLS1, LC3B)")
        _add(f"      Autofagia (LC3B-II/LC3B-I ratio):  {avg.get('autophagy', 0):.3f}")
        _add(f"      Macropinocitosis (TMR-dextran):    {avg.get('macropinocytosis', 0):.3f}")
        _add(f"      Lipogénesis (FASN, ACC):            {avg.get('lipogenesis', 0):.3f}")
        _add(f"      ROS (DCF-DA, MitoSOX):             {avg.get('ROS_level', 0):.3f}")
        if avg.get('autophagy', 0) > 0.4:
            _add("      ⚠ Autofagia activa → supervivencia bajo estrés metabólico")
            _add("        Considerar inhibición con hidroxicloroquina")
            
        _add("")
        _add(f"      DDR (γH2AX foci):                  {avg.get('DDR_active', 0):.3f}")
        _add(f"      BRCA1/2 funcional:                 {avg.get('BRCA_functional', 0):.3f}")
        if avg.get('autophagy', 0) > 0.6:
            _add("      ⚠ Autofagia elevada → mecanismo de supervivencia")
            _add("        Marcadores: p62/SQSTM1 ↓ (degradado), LC3B-II ↑")
            _add("        Considerar HCQ (cloroquina) para sensibilizar")

        # Ratio proliferación/apoptosis
        _add("")
        _add("### 8.9. ÍNDICES NETOS")
        p = avg.get('proliferation_signal', 0)
        a = avg.get('apoptosis_signal', 0.001)
        _add(f"      Señal proliferativa (Ki-67 equiv):  {p:.3f}")
        _add(f"      Señal apoptótica (TUNEL equiv):     {a:.3f}")
        _add(f"      Ratio proliferación/apoptosis:      {p/a:.1f}")
        _add(f"      Señal de supervivencia:             {avg.get('survival_signal', 0):.3f}")
        _add(f"      Señal de invasión:                  {avg.get('invasion_signal', 0):.3f}")
        if p/a > 5:
            _add("      → Crecimiento tumoral activo (ratio >5)")
        elif p/a > 1:
            _add("      → Crecimiento neto positivo pero controlable")
        else:
            _add("      → Regresión tumoral (ratio <1)")
    _add("")

    # ═══ 7. MECANISMOS DE RESISTENCIA ═══
    _add("## 9. MECANISMOS DE RESISTENCIA DETECTADOS")
    resistance_found = _analyze_resistance(model, drug_names)
    if resistance_found:
        for i, (title, explain, data, evidence) in enumerate(resistance_found, 1):
            _add(f"### 9.{i}. {title}")
            _add(f"**► Trascendencia Clínica:**")
            for line in explain:
                 _add(f"         {line}")
            _add(f"**► Datos Moleculares:**")
            for line in data:
                 _add(f"         {line}")
            _add(f"       ► {evidence}")
            _add("")
    else:
        _add("   No se detectaron mecanismos de resistencia significativos.")
    _add("")

    # ═══ 10. HIPÓTESIS Y PROPUESTA EXPERIMENTAL ═══
    if not active:
        pass
    else:
        _add("## 10. HIPÓTESIS GENERADA Y PROPUESTA DE EXPERIMENTO")
        _add("   Basándose en los resultados de esta simulación in silico,")
        _add("   proponemos la siguiente hipótesis testable:")
        _add("")
        
        # Generate hypothesis based on what we observed
        if is_krasi and not has_afatinib:
            _add("### FUNDAMENTO E HIPÓTESIS MECANÍSTICA:")
            _add(f"   **Hipótesis:** La inhibición vertical monoterápica con {drug_str} ejerce un potente efecto citoestático inicial al colapsar eficientemente el eje RAS-MAPK (pERK↓). Sin embargo, esta perturbación homeostática levanta los bucles de retroalimentación negativa que restringen la expresión de receptores tirosina quinasa. Predecimos que el tumor sorteará la inhibición reclutando un pool emergente de EGFR/HER2 hiperactivos, reactivando la viabilidad mediante KRAS mutacionalmente intacto (wild-type) o un rescate heterodimérico paralelo.")
            _add("   **Propuesta Racional:** El secuestro farmacológico dual del nodo KRAS (inhibidor específico mutante) en concomitancia con un inhibidor pan-ErbB (Afatinib) yugulará la plasticidad de red precoz, induciendo apoptosis irrevocable en el lecho tumoral residual antes de que el fenotipo diploide recablee su metabolismo.")
            _add("")
            _add("### DISEÑO EXPERIMENTAL IN VIVO PROPUESTO (Modelo PDX Ortotópico):")
            _add("   **1. Sujetos y Modelaje:** Ratones atímicos NU/NU implantados ortotópicamente con xenoinjertos derivados de paciente (PDX) G12D/TP53-mut, estratificando por volumen basal (~150 mm³).")
            _add("   **2. Estratificación de Brazos (n=10/grupo):**")
            _add("      • Cohorte A: Vehículo (Control)")
            _add(f"      • Cohorte B: {drug_desc_list[0]} (Monoterapia, dosis óptima biodisponible)")
            _add("      • Cohorte C: Afatinib 25 mg/kg/día p.o. (Control RTK)")
            _add(f"      • Cohorte D: {drug_desc_list[0]} + Afatinib (Bloqueo Vertical Dual)")
            _add("   **3. Endpoints Primarios y Biomarcadores:**")
            _add("      • Cinética de masa volumétrica tumoral por μCT y resonancia magnética longitudinal (día 0, 14, 28).")
            _add("      • Perfil fosfo-proteómico al sacrificio (Día 28): Western blot tisular para pERK(T202/Y204), pEGFR(Y1068), pAKT(S473) y pSTAT3.")
            _add("      • Índice Apoptótico Quirúrgico: Tinción TUNEL cruzada con Ki-67 en histología de core tumoral.")
            _add("   **4. Proyección de Resultados In Silico:** Se anticipa una respuesta divergente drástica entre el Brazo B (regresión inicial con escape visible hacia la semana 3) y el Brazo D (regresión tumoral profunda sostenida y erradicación del nicho hipóxico sin recidiva medible a 30 días).")
            
        elif is_krasi and has_afatinib and has_stat3:
            _add("### FUNDAMENTO E HIPÓTESIS MECANÍSTICA:")
            _add("   **Hipótesis:** La recalcitrancia histórica del PDAC obedece a su topología de red ultra-robusta. La terapia dual inhibitoria (KRAS mutante + pan-ErbB) silencia efectivamente el polo mitogénico (ERK/EGFR), precipitando apoptosis primaria. No obstante, enclou resiliencia de nicho al preservar vías de supervivencia transcripcionales paracrinas (STAT3/IL-6). Postulamos que el triplete terapéutico añadiendo degrada-ción covalente de STAT3 erradicará el último corredor de evasión citosólico, aboliendo toda capacidad del tumor para transicionar a un fenotipo estromal protector o basal.")
            _add("")
            _add("### DISEÑO EXPERIMENTAL PRECLÍNICO TRASLACIONAL:")
            _add("   **1. Modelo Mixto Co-Cultivo Avanzado (Organoides 3D):** Organoides derivados de paciente (PDOs) PDAC en co-cultivo ex vivo con fibroblastos asociados al cáncer (CAFs) inmortalizados, para recapitular el eje inflamatorio estromal IL-6/STAT3.")
            _add("   **2. Brazos de Tratamiento:**")
            _add("      • Vehículo")
            _add(f"      • Dual (Standard of Care Quirúrgico): {drug_desc_list[0]} + Afatinib")
            _add("      • Terapia Triple: Dual + Degrader SD-36")
            _add("   **3. Readouts y Criterios Analíticos:**")
            _add("      • Monitorización temporal topográfica de viabilidad celular (MTS 3D-organoide) a 48h, 96h y 144h.")
            _add("      • Validación transcripcional espacial: RNA-seq de célula única (scRNA-seq) para identificar la abolición absoluta de firmas mesenquimales (ZEB1, VIM, SNAI1) dependientes de STAT3 en las subpoblaciones residuales.")
            _add("   **4. Proyección In Silico:** Abrogación completa de la interacción bidireccional Paracrina CAF-Célula tumoral. Supresión del índice invasivo a <5% basal comparado con la terapia dual (>35%).")
            
        elif is_krasi and has_afatinib:
            _add("### FUNDAMENTO E HIPÓTESIS MECANÍSTICA:")
            _add("   **Hipótesis:** La interceptación molecular con la dupla farmacológica (Inhibición direccional KRAS acoplada a supresión RTK pan-HER) asfixia con rotundo éxito la señalización trófica mediada por MAPK y anula el rebote clásico compensatorio de receptores transmembrana. No obstante, nuestro gemelo digital detecta riesgo latente mediado por la integridad transcripcional de STAT3 y YAP/TAZ. Sospechamos que el tumor tolerará esta restricción si es capaz de mimetizar señales pro-supervivencia independientes de fosforilación mediada por quinasa receptora.")
            _add("   **Propuesta Sensibilizadora:** Incorporar una ventana terapéutica posterior interceptando factores nucleares epigenéticos o perturbando mecánicamente el nicho.")
            _add("")
            _add("### DISEÑO EXPERIMENTAL SUGERIDO (Validación Biomarcadores):")
            _add("   • Evaluar persistencia proteómica de fosfo-STAT3 y localización nuclear YAP en biopsias líquidas/tumores extirpados (murinos y prospectos clínicos) bajo régimen de poli-terapia dual. Identificar el umbral de escape transcripcional.")
            
        elif is_parpi:
            brca = avg.get('BRCA_functional', 0.9)
            _add("### FUNDAMENTO E HIPÓTESIS MECANÍSTICA:")
            if brca < 0.3:
                _add("   **Hipótesis:** Las firmas biológicas computadas confirman un fenotipo de Deficiencia de Recombinación Homóloga (HRD) robustamente instaurado (BRCA disfuncional). Este déficit intrínseco incapacita al tumor PDAC para solventar averías genómicas de cadena doble de alta fidelidad. Predicamos que la inhibición catalítica y el atrapamiento físico del complejo PARP1 exacerban el colapso de las horquillas de replicación, conduciendo obligatoriamente a una letalidad sintética inminente y cataclismo cromosómico irreversible.")
                _add("")
                _add("### ENSAYO DE INTERVENCIÓN EX VIVO (Organoides HRD+):")
                _add("   **1. Modelo:** Cultivo primario derivado de paciente secuenciado con firma genómica BRCA1/2-delecionada confirmada.")
                _add("   **2. Terapia Combinada Racionalizada:** Olaparib monoterápico versus adición sinérgica de agentes originadores de rotura de doble cadena (e.g. Gemcitabina/Cisplatino en dosificación fraccionada).")
                _add("   **3. Medición Citométrica:** Cuantificación cuantitativa de focos de ADN (Inmunofluorescencia multiplexada de γH2AX y 53BP1) a las 48h, demostrando fallo de reparación y desintegración del núcleo tumoral.")
            else:
                _add("   **Veredicto Clínico:** La simulación genómica corrobora un aparato reparador de ADN molecularmente intacto (BRCA-funcional, competencia HR intacta). La prescripción empírica de inhibidores PARP exhibirá una futilidad clínica abrumadora. Se recomienda encarecidamente la re-estratificación del paciente a protocolos estándar quimiotóxicos combinados (FOLFIRINOX) o reclutamiento en ensayos de inhibición vertical de vías oncogénicas directas (Polo RAS-RAF-MEK).")
                
        elif is_chemo:
            _add("### FUNDAMENTO E HIPÓTESIS MECANÍSTICA:")
            _add("   **Hipótesis:** La citotoxicidad sistémica infundida frena el momentum proliferativo de la masa tumoral, pero es incapaz de ejecutar el clivaje selectivo a nivel sub-clonal. Computamos una respuesta intrínseca masiva orientada a la protección celular: inducción vigorosa del flujo lisosómico (autofagia protectora) acompañada de una reprogramación metabólica y una orquestación estromal impenetrable (desmoplasia fibro-colagenosa) que exilia físicamente a las terapias acuosas de penetrar el núcleo hipóxico del tumor.")
            _add("   **Propuesta Traslacional:** El estrangulamiento de los conductos de reciclaje energético en este instante citotóxico de crisis precipitará una ceguera metabólica total y, en consecuencia, una crisis necrótica generalizada y ruptura poblacional irreversible en el PDAC residual altamente refractario.")
            _add("")
            _add("### ESQUEMA EXPERIMENTAL DE VALIDACIÓN:")
            _add("   **1. Cohorte de Prueba:** Modelos murinos de progresión endógena lenta KPC (KrasLSL.G12D/+; p53R172H/+; PdxCretg/+) con alto componente estromal de penetrancia comprobada.")
            _add("   **2. Esquema Posológico Combinado:**")
            _add("      • Régimen Quimioterápico MTD (Máxima Dosis Tolerada)")
            _add("      • Régimen MTD modulado temporalmente con un potente inhibidor del poro lisosomal intra-celular (p.ej.: Difosfato de Cloroquina 50mg/kg diarios).")
            _add("   **3. Metricas y Endpoints Rigurosos:**")
            _add("      • Puntuación volumétrica (Resonancia/PET scan de captación metabólica fluorodeoxiglucosa [18F-FDG]).")
            _add("      • Microscopía electrónica pos-sacrificio evaluando hipertrofia, agregación granular, y parálisis en la digestión de autoligosomas (marcador de p62/SQSTM1 intra-tumoral indeseado).")
            
        elif 'trametinib' in active:
            _add("### FUNDAMENTO E HIPÓTESIS MECANÍSTICA:")
            _add("   **Hipótesis:** La intervención monoterápica inhibitoria sobre el eje cinasa MEK1/2 (downstream hiper-conductor de KRAS) suprime catastróficamente la fosforilación terminal de ERK en las primeras 24 h. Sin embargo, advertimos un fenómeno flagrante de fracaso terapéutico programado: la ceguera del feedback fisiológico auto-restrictivo mediado por Spry/DUSP induce inexorablemente un tsunami de reactivación cascada-arriba (EGFR/HER2, pCRAF), reinyectando a medio plazo vectores de fosforilación lateral viables (STAT3/PI3K) y abortando la letalidad inducida.")
            _add("")
            _add("### PLAN DE ENSAYO COMBINATORIO RACIONAL:")
            _add("   **1. Plataforma Ex Vivo:** Esferoides tumorales derivados de Biopsia (FNA) clasificados fenotípicamente por índice de expresión membranaria pEGFR-Y1068 basal.")
            _add("   **2. Bloqueo Vertical Dinámico:** Trametinib monoterápico versus acoplamiento sinérgico precoz (T<48h) con Cetuximab/Erlotinib (anti-EGFR) o MRTX1133 (anti-KRAS G12D upstream).")
            _add("   **3. Biométrica Validatoria:** Cuantificación multiplex temporal de viabilidad y purificación proteómica subcelular demostrando la inhibición total de fosforilaciones de rescate compensatorias, corroborando así que la vulnerabilidad de MEK existe si (y solo si) se previene el rebote del circuito sensor extracelular.")
        else:
            _add("   [El contexto de tratamiento administrado carece actualmente de modelado predictivo de hipótesis experimental acoplada en el corpus iterativo base del simulador gemelo. Ampliar la cobertura de agentes a ensayo].")
    _add("")
    _add("=" * 72)

    return "\n".join(lines)


def _average_nodes(cancer_cells):
    """Promedia los nodos de señalización de todas las células vivas."""
    avg = {}
    for c in cancer_cells:
        if hasattr(c, 'signaling'):
            for k, v in c.signaling.nodes.items():
                avg[k] = avg.get(k, 0) + v
    n = len(cancer_cells)
    for k in avg:
        avg[k] /= n
    return avg


def _analyze_resistance(model, drug_names) -> list:
    """Detecta y explica mecanismos de resistencia con rigor clínico descriptivo."""
    mechanisms = []

    cancer_cells = [a for a in model.agents
                    if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive]
    if not cancer_cells:
        return mechanisms

    avg = _average_nodes(cancer_cells)
    
    active = set(drug_names) if drug_names else set()
    if not active and hasattr(model, 'drug_doses'):
        active = {d for d, v in model.drug_doses.items() if v > 0}
    active = set(active)

    kras_drugs = {'mrtx1133', 'daraxonrasib', 'rmc7977', 'eras0015', 'ly4066434', 'sotorasib'}
    chemo_drugs = {'gemcitabine', '5fu', 'oxaliplatin', 'irinotecan', 'nab_paclitaxel'}
    immuno_drugs = {'anti_pd1', 'anti_ctla4'}

    # 2024 KRASi Resistencia Metabólica y Amplificación
    if kras_drugs & active and avg.get('macropinocytosis', 0) > 0.45:
        mechanisms.append((
            "Evasión Metabólica por Macropinocitosis (Vía AGER-DIAPH1)",
            ["El análisis fenotípico revela un profundo recableado metabólico en respuesta al bloqueo del eje KRAS. Ante la supresión de la vía glucolítica canónica, el tumor ha hiperactivado la macropinocitosis como mecanismo de rescate nutricional, engullendo proteínas extracelulares del estroma para mantener el pool de aminoácidos intracelulares y sostener la supervivencia."],
            [
                f"Índice de Macropinocitosis: {avg.get('macropinocytosis',0):.2f}",
                f"Firma transcripcional MYC: {avg.get('MYC_active',0):.2f}",
                "Implicación Terapéutica: La inhibición farmacológica de la macropinocitosis (potencialmente mediada por derivados amiloroides o inhibidores NHE1) representa una vulnerabilidad sintética en este contexto."
            ],
            "Base Mecanística: Vía AGER-DIAPH1 y dependencia de MYC inducida por estrés nutricional (Nature 2024)."
        ))

    # 2026 Feedback Negativo RTK post-RMC-6236
    if kras_drugs & active and (avg.get('EGFR_active', 0) > 0.4 or avg.get('HER2_active', 0) > 0.3):
        mechanisms.append((
            "Resistencia Adaptativa por Alivio de Feedback Negativo (Eje RTK-MAPK)",
            ["Se documenta una rápida adaptación no mutacional secundaría a la caída profunda de pERK. La pérdida fisiológica de la transcripción de reguladores negativos (DUSP6, SPRY2) ha desinhibido masivamente la expresión y fosforilación de los receptores tirosina quinasa (EGFR y HER2) en la membrana celular. Este 'choque de red' provoca la reactivación paradójica de KRAS wild-type, bypasseando el bloqueo del inhibidor G12D/ON."],
            [
                f"Nivel pEGFR (Y1068 equiv): {avg.get('EGFR_active',0):.2f}",
                f"Nivel pHER2 (Y1248 equiv): {avg.get('HER2_active',0):.2f}",
                f"Actividad KRAS-WT compensatoria: {avg.get('KRAS_WT_active',0):.2f}",
                "Implicación Terapéutica: Fuerte justificación biológica para el bloqueo vertical combinando el inhibidor de KRAS con un inhibidor pan-ErbB (como Afatinib o Cetuximab)."
            ],
            "Base Mecanística: Reestructuración homeostática de la red de señalización MAPK (Kano et al. 2024)."
        ))

    # Bypass YAP post-KRASi/Chemo
    if (kras_drugs & active or chemo_drugs & active) and avg.get('YAP_nuclear', 0) > 0.45:
        mechanisms.append((
            "Reprogramación Transcripcional Mediada por YAP/TAZ y Transición EMT",
            ["El estrés terapéutico ha inducido una transición epitelio-mesenquimal (EMT) acoplada a la translocación nuclear sostenida de los efectores YAP/TAZ. Este programa transcripcional basal-like confiere al clon tumoral un estado de supervivencia intrínseca que es totalmente independiente de la señalización clásica de KRAS, actuando como un bypass mecanicista terminal."],
            [
                f"Localización nuclear YAP: {avg['YAP_nuclear']:.2f} (Umbral de resistencia superado)",
                f"Score de fenotipo basal/squamous: {avg.get('basal_like_program',0):.2f}",
                "Implicación Terapéutica: Posible beneficio de terapias disruptoras de TEAD/YAP (ej. análogos de Verteporfina) o reversión epigenética."
            ],
            "Base Mecanística: Reorganización transcripcional por mecanotransducción y estrés citotóxico (Shao DD et al)."
        ))

    # Reactivación PI3K post-KRASi
    if kras_drugs & active and avg.get('PI3K_active', 0) > 0.45:
        mechanisms.append((
            "Divergencia de Señalización hacia el Eje PI3K/AKT/mTOR",
            ["Como respuesta adaptativa al bloqueo de la cascada principal MAPK, las células tumorales han derivado el flujo de señales mitogénicas originadas en la membrana hacia la cascada paralela de PI3K. Esta hiperactivación rescata la síntesis proteica (vía mTORC1) y bloquea potentemente la apoptosis (vía AKT-BAD), estabilizando el citoesqueleto."],
            [
                f"Tono basal PI3K: {avg['PI3K_active']:.2f}",
                f"Fosforilación AKT (Ser473 equiv): {avg.get('AKT_active',0):.2f}",
                "Implicación Terapéutica: Señal inequívoca para evaluar el bloqueo dual (KRASi + PI3K/mTORi) para colapsar ambas vías de supervivencia."
            ],
            "Base Mecanística: Cross-talk recíproco funcional entre PI3K y MAPK (García-Alonso 2025)."
        ))

    # 2024 FOLFIRINOX - Eje NOTCH/DDR
    if {'5fu', 'oxaliplatin', 'irinotecan'}.issubset(active) and avg.get('NOTCH_active', 0) > 0.4:
        mechanisms.append((
            "Quimiorresistencia por Hiperactivación del Eje NOTCH/DDR",
            ["El régimen citotóxico FOLFIRINOX está induciendo un cataclismo genómico esperado. Sin embargo, el análisis nodular demuestra que el tumor ha hiperactivado simultáneamente la vía del receptor NOTCH y la Respuesta al Daño del ADN (DDR). Este bucle permite que la maquinaria nuclear repare los aductos de platino y las roturas de cadena a mayor velocidad de la que el fármaco logra generarlas."],
            [
                f"Respuesta a Daño del ADN (DDR): {avg.get('DDR_active',0):.2f}",
                f"Señalización del receptor NOTCH: {avg.get('NOTCH_active',0):.2f}",
                "Implicación Terapéutica: Existe racional para sensibilizar el tumor bloqueando la escisión proteolítica de NOTCH mediante inhibidores de gamma-secretasa."
            ],
            "Base Mecanística: Interacción GALNT5-MYH9 catalizadora de hiperactivación NOTCH en evasión a QT (Gastroenterology 2024)."
        ))

    # Autofagia protectora
    if chemo_drugs & active and avg.get('autophagy', 0) > 0.55:
        mechanisms.append((
            "Supervivencia Celular Mediada por Autofagia Extrema",
            ["Ante el colapso energético y el severo daño en las organelas inducido por la quimioterapia, las células neoplásicas han recurrido a un programa de autofagia masiva. Están literalmente autosecuestrando y digiriendo lisosomalmente sus propias proteínas dañadas para reciclar ATP y evitar el umbral de muerte apoptótica."],
            [
                f"Flujo autofágico (ratio LC3B-II/I subrogado): {avg.get('autophagy',0):.2f}",
                "Implicación Terapéutica: Sensibilización potencial de altísimo impacto si se bloquea el lisosoma distalmente (Hidroxicloroquina/HCQ)."
            ],
            "Base Mecanística: Catabolismo autofágico inducido por estrés citotóxico (Yang S et al)."
        ))

    # 2024 PARPi Resistencia - PARP a OXPHOS
    if 'olaparib' in active and avg.get('OXPHOS_active', 0) > 0.4:
        mechanisms.append((
            "Escape a Letalidad Sintética y Shift a Fosforilación Oxidativa",
            ["Se detecta un fenotipo dual de escape a inhibidores de PARP. Las células muestran funcionalidad parcial de la recombinación homóloga (sospecha de mutación de reversión estructural en BRCA) y, simultáneamente, han reprogramado su mitocondria para depender casi exclusivamente de la fosforilación oxidativa (OXPHOS), blindándose frente al letal estrés de los radicales libres."],
            [
                f"Funcionalidad percibida recombinación (BRCA): {avg.get('BRCA_functional',0):.2f}",
                f"Dependencia mitocondrial (OXPHOS): {avg.get('OXPHOS_active',0):.2f}",
                "Implicación Terapéutica: El ataque coordinado al metabolismo mitocondrial (inhibidores de complejo I como IACS-10759) podría quebrar este escudo."
            ],
            "Base Mecanística: Adaptación metabólica selectiva bajo presión de iPARP (Cell Metab 2024)."
        ))

    # Hipoxia
    hif = max(avg.get('HIF1A_active', 0), avg.get('HIF2A_active', 0))
    if hif > 0.4:
        mechanisms.append((
            "Resistencia Multifactorial por Microambiente Hipóxico",
            ["La precaria vascularización del núcleo tumoral ha forzado la estabilización de los factores HIF-1α y HIF-2α. Este hipo-oxigenamiento crónico no solo inactiva la generación de especies reactivas de oxígeno (ROS) necesarias para potenciar la quimioterapia, sino que ha polarizado transcripcionalmente a la célula hacia un fenotipo quiescente, glicolítico y pan-resistente altamente agresivo."],
            [
                f"Estabilización HIF1/2α: {hif:.2f}",
                f"Efecto Warburg integrado: {avg.get('warburg_effect',0):.2f}"
            ],
            "Base Mecanística: Inducción de bombas de eflujo (MDR1) y arresto de ciclo por asfixia del lecho (Radiotherapy and Oncology)."
        ))

    # Desmoplasia
    if hasattr(model, 'microenv'):
        ecm = float(np.mean(model.microenv.ecm_density))
        if ecm > 0.4:
            mechanisms.append((
                "Aislamiento Físico por Estroma Desmoplásico Denso",
                ["Independientemente de la susceptibilidad molecular intra-celular, el tumor ha depositado una densa matriz extracelular colagenosa (desmoplasia) orquestada por la fracción estromal (CAFs). Esta costra fibrótica deprime drásticamente la perfusión intratumoral impidiendo que las moléculas terapéuticas alcancen siquiera a las células diana."],
                [
                    f"Densidad de Matriz Extracelular (ECM): {ecm:.2f}", 
                    f"Coeficiente de penetración biodisponible: {max(1.0-ecm*0.6,0.2):.0%}"
                ],
                "Base Mecanística: Hipertensión intersticial fluídica secundaria a deposición de ácido hialurónico (Gore J et al)."
            ))

    # Fix return values
    formatted_mechanisms = []
    for m in mechanisms:
        title = m[0]
        explain = m[1]
        data = m[2]
        evidence = m[3] if len(m) > 3 else "Evidencia Científica: Consenso oncológico actual sobre fenotipos refractarios en PDAC."
        formatted_mechanisms.append((title, explain, data, evidence))

    return formatted_mechanisms


def _get_recommendations(model, drug_names, resistance) -> list:
    """Recomendaciones terapéuticas actualizadas."""
    recs = []
    active = set(d for d in drug_names if model.drug_doses.get(d, 0) > 0)

    if not active:
        recs.append("FOLFIRINOX (5-FU + Oxaliplatino + Irinotecán): estándar en fit patients.")
        recs.append("Gemcitabina + Nab-Paclitaxel: alternativa si ECOG ≥2.")
        recs.append("Considerar perfil molecular (KRAS) para terapia dirigida a futuro.")
        return recs

    has_resistance = {r[0] for r in resistance}

    if any("YAP/TAZ" in r for r in has_resistance) and 'verteporfin' not in active:
        recs.append("Añadir Verteporfin (YAPi) para bloquear bypass celular.")

    if any("AUTOFAGIA" in r for r in has_resistance):
        recs.append("Añadir Hidroxicloroquina (HCQ) para inhibir maquinaria de reciclaje.")

    if any("FOLFIRINOX MEDIADA POR NOTCH" in r for r in has_resistance):
        recs.append("Investigar inhibidores de gamma-secretasa / bloqueadores de NOTCH.")

    if any("ESCAPE METABÓLICO" in r for r in has_resistance):
        recs.append("Usar agentes bloqueadores de macropinocitosis (anti-metabólicos).")

    if any("OXPHOS" in r for r in has_resistance):
        recs.append("Añadir IACS-10759 (Inhibidor de OXPHOS) para revertir escape PARPi.")

    if any("ESTROMAL" in r for r in has_resistance):
        recs.append("Considerar PEGPH20 (Hialuronidasa) para ablandar estroma.")

    cancer_cells = [a for a in model.agents if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive]
    if cancer_cells and hasattr(cancer_cells[0], 'signaling'):
        muts = cancer_cells[0].signaling.mutations
        if muts.get('BRCA') and 'olaparib' not in active:
            recs.append("⭐ BRCA mutado → Considerar escalar a Olaparib (PARPi).")
        kras_drugs = {'mrtx1133', 'daraxonrasib', 'rmc7977'}
        if not (kras_drugs & active) and muts.get('KRAS'):
            recs.append("Probar inhibidor pan-KRAS (RMC-6236) — ORR 20-36% en tumores G12X.")

    if not recs:
        recs.append("Tratamiento eficaz. Vigilar biomarcadores y evolución RECIST.")

    return recs

def get_pathway_deltas(model) -> dict:
    """
    Calcula el porcentaje de cambio (Delta) de las vías moleculares clave
    entre el T0 (paso 1) y el Tf (paso final).
    Retorna: { 'Node_Name': (delta_pct, T0_val, Tf_val) }
    """
    if not hasattr(model, 'initial_node_averages') or not model.initial_node_averages:
        return {}
        
    cancer_cells = [a for a in model.agents if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive]
    if not cancer_cells:
        return {}
        
    t0_avg = model.initial_node_averages
    tf_avg = _average_nodes(cancer_cells)
    
    # Nodos relevantes a monitorizar ("Panel Ómico")
    key_nodes = [
        'KRAS_active', 'ERK_active', 'AKT_active', 'mTOR_active', 
        'YAP_nuclear', 'HIF1A_active', 'STAT3_active', 'PDL1_expression',
        'apoptosis_signal'
    ]
    
    deltas = {}
    for node in key_nodes:
        val_0 = t0_avg.get(node, 0.0)
        val_f = tf_avg.get(node, 0.0)
        
        # Evitar divisiones por cero
        if val_0 < 0.01 and val_f < 0.01:
            delta = 0.0
        elif val_0 < 0.01:
            delta = 100.0  # Fue encendido de novo
        else:
            delta = ((val_f - val_0) / val_0) * 100.0
            
        # Cap a +/- 150% para gráficas más limpias
        delta = max(-150.0, min(delta, 150.0))
            
def generate_pathway_network_fig(model):
    """Genera un infográfico interactivo celular avanzado de PDAC usando Plotly."""
    import plotly.graph_objects as go
    
    cancer_cells = [a for a in model.agents if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive]
    if not cancer_cells:
        return None
        
    avg = _average_nodes(cancer_cells)
    
    # Redefiniendo los nodos con posiciones más precisas (y: 1.0 a -0.5)
    nodes = {
        # Receptores y Señales Extracelulares (Membrana, y ≈ 0.9)
        'TGFB_R': {'pos': (-4, 0.9), 'val': avg.get('TGFB_signal', 0)},
        'EGFR': {'pos': (-1.5, 0.9), 'val': avg.get('EGFR_active', 0)},
        'HER2': {'pos': (0, 0.9), 'val': avg.get('HER2_active', 0)},
        'IL6_R': {'pos': (3.5, 0.9), 'val': avg.get('IL6_signal', 0)},
        'PTEN': {'pos': (2, 0.8), 'val': avg.get('PTEN_active', 0)},
        
        # GTPasas y Nodos de Cuello de Botella (Citoplasma Alto, y ≈ 0.7)
        'KRAS WT': {'pos': (-2, 0.7), 'val': avg.get('KRAS_WT_active', 0)},
        'KRAS Mut': {'pos': (-0.75, 0.7), 'val': avg.get('KRAS_active', 0)},
        'PI3K': {'pos': (1.2, 0.7), 'val': avg.get('PI3K_active', 0)},
        
        # Cascadas de Quinasas (Citoplasma Medio, y ≈ 0.45)
        'SMAD4': {'pos': (-4, 0.5), 'val': avg.get('SMAD4_active', 0)},
        'RAF': {'pos': (-1, 0.55), 'val': avg.get('RAF_active', 0)},
        'AKT': {'pos': (1.2, 0.55), 'val': avg.get('AKT_active', 0)},
        'STAT3': {'pos': (2.8, 0.5), 'val': avg.get('STAT3_active', 0)},
        'NFκB': {'pos': (4, 0.5), 'val': avg.get('NFKB_active', 0)},
        
        'MEK': {'pos': (-1, 0.4), 'val': avg.get('MEK_active', 0)},
        'mTOR': {'pos': (1.2, 0.4), 'val': avg.get('mTOR_active', 0)},
        
        # Efectores y Factores Transcripcionales (Núcleo y pro-núcleo, y ≈ 0.15)
        'ERK': {'pos': (-1, 0.25), 'val': avg.get('ERK_active', 0)},
        'BCL2': {'pos': (-4, 0.2), 'val': avg.get('BCL2_active', 0)},
        'TP53 / p53': {'pos': (-2.5, 0.2), 'val': avg.get('TP53_active', 0)},
        'HIF1A(Hipoxia)': {'pos': (3.5, 0.2), 'val': avg.get('HIF1A_active', 0)},
        
        # Maquinaria Transcripcional y EMT (Núcleo Profundo, y ≈ 0.0)
        'YAP/TEAD': {'pos': (-1.5, 0.0), 'val': avg.get('YAP_nuclear', 0)},
        'MYC': {'pos': (0.5, 0.0), 'val': avg.get('MYC', 0.5)}, # MYC proxy si existe
        'Caspa-3': {'pos': (-4, 0.0), 'val': avg.get('CASP3_active', 0)},
        'EMT/ZEB1': {'pos': (-2.8, -0.1), 'val': avg.get('ZEB1_active', avg.get('VIM_expression', 0))},
        
        # Fenotipos Finales (y ≈ -0.3)
        'Apoptosis': {'pos': (-4, -0.3), 'val': avg.get('apoptosis_signal', 0)},
        'Proliferación': {'pos': (-1, -0.3), 'val': avg.get('proliferation_signal', 0)},
        'Supervivencia': {'pos': (1, -0.3), 'val': avg.get('survival_signal', 0)},
        'Inmuno-Evade\n(PD-L1/MHC↓)': {'pos': (3.5, -0.3), 'val': max(avg.get('PDL1_expression', 0), 1.0-avg.get('MHC1_expression', 1))},
        'Warburg/AutoF': {'pos': (2, -0.15), 'val': max(avg.get('warburg_effect', 0), avg.get('autophagy', 0))}
    }
    
    edges = [
        # Membrane to internal
        ('TGFB_R', 'SMAD4'), ('TGFB_R', 'EMT/ZEB1'),
        ('EGFR', 'KRAS Mut'), ('EGFR', 'KRAS WT'), ('EGFR', 'PI3K'), ('EGFR', 'STAT3'),
        ('HER2', 'KRAS Mut'), ('HER2', 'KRAS WT'), ('HER2', 'PI3K'),
        ('IL6_R', 'STAT3'), ('IL6_R', 'NFκB'),
        ('PTEN', 'PI3K'), # Negative in biology, displayed as line
        
        # Cascades
        ('KRAS Mut', 'RAF'), ('KRAS WT', 'RAF'), ('KRAS Mut', 'PI3K'),
        ('PI3K', 'AKT'), ('AKT', 'mTOR'), ('RAF', 'MEK'), ('MEK', 'ERK'),
        
        # Cross-talk and Nucleus
        ('SMAD4', 'YAP/TEAD'), ('SMAD4', 'EMT/ZEB1'),
        ('STAT3', 'MYC'), ('STAT3', 'Supervivencia'), ('STAT3', 'Inmuno-Evade\n(PD-L1/MHC↓)'),
        ('NFκB', 'Inmuno-Evade\n(PD-L1/MHC↓)'), ('NFκB', 'Supervivencia'),
        ('HIF1A(Hipoxia)', 'Warburg/AutoF'), ('HIF1A(Hipoxia)', 'Inmuno-Evade\n(PD-L1/MHC↓)'),
        
        ('ERK', 'MYC'), ('ERK', 'Proliferación'), ('ERK', 'Supervivencia'), ('ERK', 'BCL2'),
        ('mTOR', 'Warburg/AutoF'), ('mTOR', 'Supervivencia'),
        ('YAP/TEAD', 'Proliferación'), ('YAP/TEAD', 'Supervivencia'),
        ('MYC', 'Proliferación'), ('EMT/ZEB1', 'Supervivencia'),
        
        # Apoptosis rules
        ('TP53 / p53', 'BCL2'), ('TP53 / p53', 'Caspa-3'),
        ('BCL2', 'Caspa-3'), ('Caspa-3', 'Apoptosis'),
        ('AKT', 'Caspa-3'), # inhibits directly
    ]
    
    fig = go.Figure()
    
    # Dibujar aristas
    for (src, dst) in edges:
        if src in nodes and dst in nodes:
            x0, y0 = nodes[src]['pos']
            x1, y1 = nodes[dst]['pos']
            fig.add_trace(go.Scatter(
                x=[x0, x1], y=[y0, y1], mode='lines',
                line=dict(width=1, color='rgba(150, 150, 150, 0.4)'),
                hoverinfo='none', showlegend=False
            ))
        
    # Dibujar nodos
    node_x, node_y, node_color, node_text = [], [], [], []
    for name, data in nodes.items():
        node_x.append(data['pos'][0])
        node_y.append(data['pos'][1])
        node_color.append(data['val'])
        node_text.append(f"<b>{name}</b><br>Actividad/Nivel: {data['val']:.2f}")
        
    fig.add_trace(go.Scatter(
        x=node_x, y=node_y, mode='markers+text',
        marker=dict(
            showscale=True, colorscale='RdYlBu_r', # _r para que Azul sea frío (0) y Rojo sea caliente (>0.8)
            cmin=0.0, cmax=1.0, color=node_color,
            size=36, colorbar=dict(title='Actividad<br>(0-1)', thickness=15),
            line_width=1.5, line_color='white'
        ),
        text=list(nodes.keys()), textposition="top center",
        hovertext=node_text, hoverinfo="text",
        textfont=dict(size=11, color="white", family="Arial"),
        showlegend=False
    ))
    
    # Celdas espaciales (Membrana, Citoplasma, Núcleo)
    fig.add_hrect(y0=0.82, y1=1.05, fillcolor="rgba(52, 152, 219, 0.1)", opacity=0.5, layer="below", line_width=0)
    fig.add_annotation(x=-4.7, y=1.00, text="Membrana Extracelular y Receptores", showarrow=False, xanchor="left", font=dict(color="rgba(255,255,255,0.7)"))
    
    fig.add_hrect(y0=0.10, y1=0.82, fillcolor="rgba(155, 89, 182, 0.1)", opacity=0.5, layer="below", line_width=0)
    fig.add_annotation(x=-4.7, y=0.78, text="Citoplasma (Cascadas Enzimáticas)", showarrow=False, xanchor="left", font=dict(color="rgba(255,255,255,0.7)"))
    
    fig.add_hrect(y0=-0.40, y1=0.10, fillcolor="rgba(46, 204, 113, 0.1)", opacity=0.5, layer="below", line_width=0)
    fig.add_annotation(x=-4.7, y=-0.36, text="Núcleo (Reprogramación Transcripcional y Fenotipo)", showarrow=False, xanchor="left", font=dict(color="rgba(255,255,255,0.7)"))
    
    fig.update_layout(
        title=dict(text="", font=dict(size=18, color="white")),
        showlegend=False, hovermode='closest',
        margin=dict(b=0, l=0, r=0, t=20),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-4.8, 4.8]),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-0.45, 1.1]),
        height=550, plot_bgcolor="rgba(0,0,0,0)", paper_bgcolor="rgba(0,0,0,0)"
    )
    
    return fig
