# Bobo Plan2BOQ

แบบ (รูปภาพ / PDF / DXF) → วิเคราะห์ผนัง-ช่องเปิด-ห้อง → โครง คสล. เบื้องต้น (ACI 318-19) → โมเดล 3D → ถอดแบบ (QTO) → BOQ ไทย 5 หมวด

- ใช้งาน: เปิด `index.html` (ไฟล์เดียว ทำงานในเบราว์เซอร์ทั้งหมด ไม่ส่งไฟล์แบบขึ้น server)
- `src/engine.js` ส่วนคำนวณล้วน (ไม่มี DOM) · `src/ui.html` ส่วน UI (three.js r128, pdf.js 3.11)
- Build: `python3 src/build.py` (รวม engine + ตัวอย่างเข้า index.html)
- หน่วย: cm, kgf, ksc · ผลออกแบบเป็น preliminary ต้องตรวจสอบด้วยการวิเคราะห์โครงสร้างจริง

ตัวอย่างทดสอบจาก: grebtsew/FloorplanToBlender3d, omarsamy3/DXF_To_Revit
หลักการจาก: datadrivenconstruction/DDC_Skills_for_AI_Agents_in_Construction (MIT)
