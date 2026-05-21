역할: 너는 서울대학교 생물정보학 및 실습 1 과목의 텀 프로젝트(Term Project) 지도를 맡은 바이오인포매틱스 분야의 수석 연구원이자 교수님이다. 
내가 제공한 1~3주차 Guided Mission 노트북 실습 내용(featureCounts 발현량 정량, bedtools 구간 연산, 단일 염기 Shannon Entropy 계산 및 UCSC 시각화)과 제공된 데이터(RNA-seq, RPF, CLIP-seq)를 완벽하게 이해하고 있다.

[요청 사항]
이번 텀 프로젝트의 최종 발표(Presentation) 및 최종 보고서에서 연구의 완결성과 독창성을 증명하기 위해, 기존 실습 데이터들을 유기적으로 결합하여 수행할 수 있는 '최종 확장 실습 주제 3가지'를 추천해 주라. 학계 논문 수준의 분석이 되도록 각 주제별 '목적 및 생물학적 배경', '상세 실행 방법론(Step-by-step Methods)', 그리고 '예상되는 시각화(Figure) 및 발표 스토리라인'을 매우 구체적으로 작성해 주라.

---

### 🚀 확장 실습 1: Lin28a Knockdown에 따른 유전자별 번역 효율(Translation Efficiency, TE) 변동 및 세포 내 위치(Localization) 상관관계 규명
(1주차 featureCounts 정량 + 2주차 RPF 데이터 처리 + 외부 데이터 결합)

1. 목적 및 배경:
   - mRNA의 양(RNA-seq)과 실제 번역되는 양(RPF)은 일치하지 않는다. Lin28a가 억제되었을 때, 어떤 유전자 군의 번역 효율(TE)이 유의미하게 변화하는지 정량화한다.
   - 1주차에 사용한 마우스 단백질 세포 내 위치(Cellular Localization) 주석 데이터를 결합하여, 특정 위치(예: 세포질 vs 핵)를 타겟팅하는 유전자들이 Lin28a의 제어를 더 강하게 받는지 통계적으로 검증한다.
2. 상세 실행 방법론 (Methods):
   - Step 1: featureCounts를 사용해 RNA-control(siLuc), RNA-siLin28a, RPF-control(siLuc), RPF-siLin28a 4개 조건의 유전자별 raw read count 행렬 생성.
   - Step 2: 총 read 수 차이를 보정하기 위해 pandas를 사용하여 각 유전자의 카운트를 TPM(Transcripts Per Million) 또는 CPM으로 정규화(Normalization).
   - Step 3: 각 유전자별 번역 효율 계산 식 구현: TE = (RPF 정규화 값) / (RNA 정규화 값)
   - Step 4: Control 대비 siLin28a 조건에서의 TE 변화량(Fold Change) 계산: Delta TE = log2(TE_siLin28a) - log2(TE_siLuc)
   - Step 5: mouselocalization-20210507.txt 데이터를 유전자 ID 기준으로 병합하고, 단백질 위치 그룹별(Cytoplasm, Nucleus 등) Delta TE 값의 분포 차이를 ANOVA 또는 t-test로 통계적 유의성 검정.
3. 시각화 및 발표 활용 포인트:
   - X축 log2(RNA 변화량), Y축 log2(RPF 변화량)인 Scatter plot을 그리고 TE가 급격히 변한 유전자 강조.
   - 세포 내 위치별 Delta TE 분포를 보여주는 Box plot 또는 Violin plot 제시. "Lin28a는 특정 위치의 유전자 번역을 선택적으로 조절한다"는 강력한 발표 결론 도출.

---

### 🚀 확장 실습 2: 개시 코돈(Start Codon) 주변 리보솜 주기성(Ribosome Periodicity) 및 Triplet Phasing의 정량적 검증
(2주차 bedtools를 이용한 개시 코돈 주변 5' 말단 카운트 가공 고도화)

1. 목적 및 배경:
   - 리보솜은 mRNA 위를 3nt(1코돈)씩 불연속적으로 이동하며 번역틀(Reading Frame)을 유지한다. RPF 데이터의 5' 말단 분포에서 나타나는 3bp 주기의 신호(Triplet Phasing)를 데이터로 증명하여 시퀀싱 데이터의 퀄리티와 리보솜 프로파일링의 생물학적 기전을 확인한다.
2. 상세 실행 방법론 (Methods):
   - Step 1: 2주차에서 `bedtools intersect`로 필터링한 `fivepcounts-filtered-RPF-siLuc.txt` 데이터를 파이썬으로 로드.
   - Step 2: 유전체의 절대 좌표를 개시 코돈(Start Codon)의 첫 번째 염기 서열 자리를 '0'으로 하는 상대 좌표(Relative Coordinate, 예: -30bp ~ +90bp)로 매핑 변환.
   - Step 3: 모든 유전자에 대해 상대 좌표 포지션별 리드 카운트를 누적 합산(Aggregation)하여 Meta-gene count 테이블 구축.
   - Step 4: 개시 코돈 이후인 +0, +3, +6, +9... 자리와 +1, +4, +7... 및 +2, +5, +8... 자리의 리드 밀도(Frame 0, 1, 2)를 비교하여 특정 프레임에 리드가 편향되어 있는지 통계적 검정(Chi-square test).
   - Step 5: 주기성을 정밀 검증하기 위해 파이썬 `scipy.signal`을 이용하여 데이터의 Autocorrelation(자기상관분석)을 수행하거나 FFT(Fast Fourier Transform)를 적용하여 '주기 3' 신호의 파워 스펙트럼 피크 검출.
3. 시각화 및 발표 활용 포인트:
   - X축을 개시 코돈으로부터의 거리(bp), Y축을 누적 리드 카운트로 하는 Meta-gene profile 바 차트(Bar chart) 작성. 3bp 간격으로 솟구치는 피크를 시각적으로 선명하게 제시하여 리보솜 프로파일링 실험의 엄밀성을 증명.

---

### 🚀 확장 실습 3: CLIP-seq 데이터를 활용한 Lin28a 결합 부위의 Shannon Entropy 분석 및 서열 모티프(Motif) 역추적
(3주차 단일 염기 수준 변이율 계산 알고리즘 및 시각화 확장)

1. 목적 및 배경:
   - CLIP-seq 분석 시 단백질과 RNA가 자외선(UV)으로 교차결합(Cross-linking)된 지점은 역전사 과정에서 염기 변형이나 결실 오류가 집중적으로 발생한다.
   - 3주차에 구현한 Shannon Entropy 계산 방식을 CLIP-seq 데이터에 적용하여, 엔트로피 피크가 발생하는 좌표를 역추적함으로써 Lin28a가 실제로 결합하는 RNA 상의 표적 서열 모티프(Motif)를 직접 발굴한다.
2. 상세 실행 방법론 (Methods):
   - Step 1: CLIP-35L33G.bam 파일로부터 리드가 많이 매핑된 주요 타겟 유전자 구간을 스캔하여, 각 포지션별 A, C, G, T, Deletion 카운트 계산.
   - Step 2: 3주차 수식을 확장하여 단일 염기 레벨 Shannon Entropy 값을 전수 계산하고 DataFrame화.
   - Step 3: 엔트로피가 주변 백그라운드 대비 유의미하게 높은 지점(예: 상위 1% 피크 혹은 특정 Threshold 이상)의 유전체 좌표(Chr, Start, End)를 추출하여 BED 파일 형태로 저장.
   - Step 4: 추출된 변이 피크 좌표를 중심으로 좌우 15bp 영역의 레퍼런스 유전체 서열(mm39 fasta)을 `bedtools getfasta` 또는 Biopython을 활용해 FASTA 파일로 일괄 추출.
   - Step 5: 추출된 Fasta 서열들을 MEME Suite(DREME) 웹 툴에 입력하거나, 파이썬 코드로 직접 4-mer/5-mer 빈도를 계산하여 기존 논문에 알려진 Lin28a 결합 서열(예: GGAGA 패턴)이 유의미하게 농축(Enrichment)되어 있는지 검증.
3. 시각화 및 발표 활용 포인트:
   - UCSC Genome Browser에 우리가 찾은 엔트로피 피크 트랙과 유전자 모델을 매핑한 스크린샷 제시.
   - 최종적으로 발굴된 핵심 결합 서열의 WebLogo 시각화 결과를 제시하여, "컴퓨터 연산(엔트로피)만으로 실제 단백질이 결합하는 분자 생물학적 타겟 서열을 찾아냈다"는 텀 프로젝트의 하이라이트 구성.

---

### ❓ 심사위원/교수님 예상 질문 디펜스 리스트 (Q&A 10가지)
추후 텀 프로젝트 발표 평가 시, 심사위원 및 교수님께서 분석 방법론과 생물학적 의미에 대해 질문할 수 있는 날카로운 질문 10가지를 엄선하고 이에 대한 완벽한 답변(Rationale)을 추가해 주라. 

(질문 리스트 예시 프로토타입)
1. 발현량 정규화 시 RPKM/FPKM 대신 왜 TPM 혹은 CPM을 사용했는가? 그 통계적 차이는 무엇인가?
2. RNA-seq과 RPF의 리드 카운트 비율만으로 TE(번역 효율)를 정의할 때 발생할 수 있는 노이즈나 한계점은 무엇이며 이를 어떻게 보완했는가?
3. RPF 데이터에서 5' 말단 포지션만 카운트(bedtools -5)하여 주기성을 보았는데, 왜 3' 말단이나 리드의 센터(Center) 대신 5' 말단을 기준으로 잡아야 리보솜의 위치 정보를 더 잘 대변하는가?
4. Triplet Phasing 분석 시 프레임(0, 1, 2) 간 차이를 Chi-square test로만 검정해도 충분한가? Autocorrelation이나 FFT를 도입했을 때의 추가적인 장점은 무엇인가?
5. CLIP-seq에서 Shannon Entropy가 높게 나오는 현상이 실제 단백질 결합에 의한 Cross-linking 오류인지, 아니면 유전체 자체의 Polymorphism(SNP)이나 시퀀싱 자체 오류인지 어떻게 구별할 수 있는가?
6. 엔트로피 피크 주변의 서열 윈도우 크기를 하필 좌우 15bp로 설정한 근거는 무엇인가? 너무 넓거나 좁을 때 어떤 문제가 생기는가?
7. Lin28a 결합 모티프를 발굴할 때 백그라운드(Control) 서열 세트는 어떻게 구성했는가? 모티프의 유의성을 증명하기 위해 어떤 통계적 enrichment 분석을 썼는가?
8. 데이터 처리 과정에서 미미하게 발현되는(Read Count가 매우 낮은) 유전자들이 TE 계산이나 엔트로피 계산 시 왜곡을 일으킬 수 있는데, 이를 필터링하기 위해 어떤 Cut-off 기법을 적용했는가?
9. 생물학적 관점에서 Lin28a 넉다운 시 특정 세포 내 위치(Localization) 그룹의 번역 효율이 대폭 변했다면, 이것이 Lin28a 단백질 고유의 기능과 어떻게 연결되는가? (기존 문헌과의 정합성)
10. 이 3가지 확장 분석 파이프라인을 구축하면서 직면했던 가장 큰 드라이-랩(Dry-lab) 측면의 데이터 병목 현상이나 컴퓨팅 이슈는 무엇이었고 이를 코딩 측면에서 어떻게 최적화했는가?

---

위의 3가지 확장 실습 주제를 구체적으로 수행할 수 있는 파이썬(Python) 및 리눅스(Linux) 핵심 소스 코드를 제공해 주고, 위에서 제시된 10가지 예상 질문에 대해 교수님이 만족하실 수밖에 없는 논리적이고 정밀한 생물정보학적 답변을 완벽하게 디테일을 채워서 작성해 주라.

---

## 📋 과제 제출 안내
**과제 제출일은 GitHub에 올라온 commit 시간 기준이 아니라 push 시간을 기준으로 합니다.**

이 안내를 통해 모든 commit 및 push를 정확하게 관리할 수 있습니다. 최종 제출 시간은 GitHub 저장소에 마지막으로 push된 시간이 기준이 되므로, 제출 마감 시간 전에 반드시 push를 완료해야 합니다.