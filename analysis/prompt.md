analysis1 파일에 있는 3 노트북 실습의 각 목적 및 실습의 method 를 요약해서 exp_summarize.md 에 작성해줘. 이 3 실습 내용 토대로 추가로 할 수 있는 실습을 해야하니 최대한 자세히 설명해줘

업로드해주신 3가지 **Guided Mission 노트북(1주차~3주차)** 실습 내용을 바탕으로, 각 실습의 목적(Objective)과 방법론(Method)을 상세히 요약해 드립니다. 아울러 이 지식을 확장하여 텀 프로젝트나 추가 연구로 진행할 수 있는 **심화 실습 아이디어**를 제안합니다.

---

### 1. 1주차 실습 (Guided Mission 1)

#### 🎯 실습 목적

* **환경 구축**: 대용량 차세대 염기서열 분석(NGS) 데이터를 다루기 위한 리눅스 명령어 기반의 분석 환경(Bioconda)을 Google Colab에 세팅합니다.
* **발현량 정량화(Quantification)**: 매핑된 시퀀싱 데이터(`.bam`)와 유전체 주석 파일(`.gtf`)을 결합하여, 유전자별(Gene-level) Read Count 테이블을 직접 가공하고 도출합니다.
* **데이터 탐색 및 병합**: 정량화된 발현량 데이터와 외부의 단백질 세포 내 위치(Cellular Localization) 주석 데이터를 결합하여 데이터 분포를 시각화할 준비를 합니다.

#### 🛠 실습 Method

1. **Bioconda 및 도구 설치**: `colab-biolab` 저장소를 복제하고 콘다 환경을 활성화한 후, Read 카운팅 프로그램인 `subread` 패키지를 설치합니다.
2. **데이터 무결성 검증**: 다운로드한 대용량 BAM 파일들의 `md5sum` 체크섬을 비교하여 파일 손상 여부를 검증합니다.
3. **FeatureCounts 실행**: `featureCounts` 명령어를 사용하여 유전자 영역(`-t exon -g gene_id`)에 오버랩되는 시퀀싱 리드 수를 조건별(`RNA-control`, `RNA-siLin28a`, `RNA-siLuc` 등)로 카운트하여 `read-counts.txt` 행렬을 생성합니다.
4. **외부 데이터 결합**: 웹 서버에서 마우스 단백질 위치 데이터(`mouselocalization-20210507.txt`)를 `pandas`로 불러와 유전자 ID(`gene_id`)를 기준으로 발현 데이터와 병합(Merge)합니다.

---

### 2. 2주차 실습 (Guided Mission 2)

#### 🎯 실습 목적

* **리보솜 프로파일링(Ribosome Profiling, RPF) 데이터 분석**: 전사체(RNA-seq) 데이터와 번역체(RPF) 데이터를 비교하여 실제 번역(Translation)이 일어나는 양상을 추적합니다.
* **유전체 브라우저 데이터 준비**: 유전자의 개시 코돈(Start Codon) 주변에서 리보솜이 어떻게 위치하는지 분포(Ribosome Footprint Meta-gene plot)를 확인하기 위해 5' 말단 count 데이터를 필터링하고 가공합니다.

#### 🛠 실습 Method

1. **Genomic Interval 가공**: `bedtools`를 활용할 수 있도록 유전체 주석 정보에서 개시 코돈을 포함하는 엑손(Exon) 영역만 추출하여 별도의 BED 파일(`gencode-exons-containing-startcodon.bed`)로 추출합니다.
2. **5' 말단 Coverage 계산**: 리보솜이 위치한 정확한 자리를 특정하기 위해, `bedtools genomecov` 명령어의 `-5` 옵션을 주어 RPF 리드의 5' 말단 포지션만을 카운트하여 `bedgraph` 형태(`fivepcounts-RPF-siLuc.bed`)로 변환합니다.
3. **공간적 교집합 필터링**: `bedtools intersect`를 사용하여 5' 말단 카운트 데이터 중에서 오직 '개시 코돈을 포함하는 엑손 영역'과 겹치는 리드만 골라내어 연산 효율성을 높이고 타겟 분석을 수행합니다.

---

### 3. 3주차 실습 (Guided Mission 3)

#### 🎯 실습 목적

* **시퀀싱 오류 및 변이 탐색 (Shannon Entropy)**: 교차 매핑(Cross-linking)이나 RNA 구조, 혹은 분석 과정에서의 염기 변형 등으로 인해 특정 유전체 위치에서 발생하는 서열 다양성을 통계학적으로 수치화합니다.
* **UCSC 유전체 브라우저 시각화**: 계산된 Shannon Entropy 값을 트랙(Custom Track) 형태로 변환하여 실제 유전체 브라우저 상에서 시각적으로 감상하고 분석합니다.

#### 🛠 실습 Method

1. **염기별 카운팅 (Base-resolution counting)**: BAM 파일로부터 특정 타겟 유전자 포지션의 각 염기(A, C, G, T) 개수를 자리별로 정밀하게 집계합니다.
2. **정보 엔트로피 계산**: 각 자리에 존재하는 염기 비율($p_i$)을 기반으로 Shannon Entropy 식($-\sum p_i \log_2 p_i$)을 파이썬 코드로 구현하여 계산합니다. (염기가 100% 일정하면 엔트로피는 0이며, 여러 염기가 섞여 있을수록 엔트로피 값이 커집니다.)
3. **Bedgraph 포맷 출력**: 유전체 좌표(`Chr`, `Start`, `End`) 뒤에 계산한 엔트로피 점수를 붙여 4열짜리 `bedgraph` 표준 규격 파일로 생성합니다.
4. **UCSC 웹 브라우저 업로드**: 생성한 파일을 UCSC Genome Browser(mm39 마우스 빌드)의 'Add Custom Tracks' 메뉴를 통해 업로드하고 유전자 모델 위에 매핑된 엔트로피 피크를 분석합니다.

---

### 💡 이 내용을 토대로 확장할 수 있는 추가 실습 아이디어

위 3가지 실습은 유전자 발현 정량화(1주차), 유전체 구간 연산 및 메타 유전자 플롯팅(2주차), 단일 염기 수준 통계 계산 및 브라우저 시각화(3주차)의 핵심 파이프라인을 관통하고 있습니다. 이를 결합하여 다음과 같은 고도화된 실습을 설계할 수 있습니다.

#### 아이디어 A: Lin28a 넉다운(siLin28a)에 따른 번역 효율(Translation Efficiency) 변화 분석 (1주차 + 2주차 확장)

* **목적**: `Lin28a` 단백질이 억제되었을 때, 세포 내 전체 유전자의 mRNA 양(RNA-seq) 변화와 실제 단백질로 만들어지는 양(RPF)의 변화가 일치하는지 규명합니다.
* **방법**:
1. 1주차에서 배운 `featureCounts`를 이용해 `RNA-siLuc`(컨트롤), `RNA-siLin28a`, `RPF-siLuc`, `RPF-siLin28a` 데이터의 유전자별 카운트를 전부 구합니다.
2. 각 샘플의 Read Count를 총 Read 수로 보정(CPM 혹은 TPM 정규화)합니다.
3. 번역 효율(TE, Translation Efficiency)을 $\text{TE} = \frac{\text{RPF Count}}{\text{RNA Count}}$로 정의하여 각 유전자마다 계산합니다.
4. 컨트롤 상태의 TE와 Lin28a가 넉다운된 상태의 TE를 비교하여, Lin28a가 주로 어떤 유전자 군의 번역을 억제하거나 촉진하는지 스캐터 플롯(Scatter plot)을 그리고, 1주차에 했던 Localization(세포 내 위치) 데이터와 융합하여 "세포질에 있는 유전자가 번역 유도가 더 잘 되는지" 등의 상관관계를 분석합니다.



#### 아이디어 B: 리보솜 주기성(Ribosome Periodicity) 및 Triplet Phasing 분석 (2주차 확장)

* **목적**: 리보솜은 mRNA 위를 코돈 단위(3nt씩)로 이동하므로, RPF 데이터는 개시 코돈 이후 3뉴클레오타이드 주기로 리드가 몰리는 특징(Periodicity)이 있습니다. 이를 데이터로 확인합니다.
* **방법**:
1. 2주차에서 구한 개시 코돈 주변 `fivepcounts-filtered-RPF-siLuc.txt` 데이터를 가져옵니다.
2. 개시 코돈의 첫 번째 염기 서열 자리를 `0`으로 잡고, RPF 5' 말단이 위치한 상대적 거리(-50bp ~ +100bp)를 구합니다.
3. 상대적 위치별로 리드 카운트를 합산(Aggregation)하여 메타 프로롯(Meta-gene profile)을 그립니다.
4. 개시 코돈 이후에 3bp 간격으로 피크(Phasing)가 나타나는지 확인하고, 주기성을 더 명확히 보기 위해 푸리에 변환(FFT)이나 자동상관분석(Autocorrelation)을 적용해 봅니다.



#### 아이디어 C: CLIP-seq 데이터를 활용한 Lin28a의 결합 모티프(Binding Motif) 및 변이율 분석 (3주차 확장)

* **목적**: 데이터팩에 포함된 `CLIP-35L33G.bam` 데이터는 Lin28a 단백질이 RNA에 직접 결합한 부위를 캡처한 것입니다. Lin28a가 결합하는 부위의 특성을 Shannon Entropy와 서열 분석으로 알아냅니다.
* **방법**:
1. 3주차 실습 방법론을 확장하여 `CLIP-35L33G.bam` 전체 혹은 특정 chromosome 영역에 대해 염기 카운팅을 수행하고 Shannon Entropy를 계산합니다.
2. CLIP-seq 특성상 단백질과 RNA가 자외선(UV)으로 교차결합(Cross-linking)된 자리는 역전사 과정에서 오류(Mutation 혹은 Deletion)가 자주 발생하여 Shannon Entropy가 높게 측정됩니다.
3. 엔트로피 값이 유독 높은 포지션(피크 부위)의 유전체 좌표들을 추출합니다.
4. 해당 좌표 주변 ±20bp 서열을 Fasta 파일로 추출한 후, **MEME Suite** 같은 모티프 분석 툴에 입력하거나 파이썬으로 4-mer/5-mer 빈도를 계산하여 Lin28a가 선호하는 표적 서열 모티프(예: GGAGA 등)를 직접 발굴해 봅니다.

앞서 요약해 드린 1~3주차 실습(유전자 발현 정량화, 유전체 구간 연산 및 메타 유전자 분석, 단일 염기 수준 변이 및 브라우저 시각화)의 데이터와 방법론을 유기적으로 연결하여, 학기 말 텀 프로젝트(Term Project)나 최종 보고서로 발전시키기 좋은 **3가지 확장 실습 주제**를 추천해 드립니다.

각 주제는 기존 실습에서 사용한 데이터(`RNA-seq`, `RPF`, `CLIP-seq`)의 특성을 극대화하고 새로운 생물학적 인사이트를 도출할 수 있도록 구성했습니다.

---

### 🚀 추천 1: Lin28a 넉다운(siLin28a)에 따른 유전자별 번역 효율(Translation Efficiency, TE) 변화 분석

* **연계 실습**: 1주차(featureCounts 발현량 정량) + 2주차(RPF 데이터 처리)

#### 1) 실습 목적

세포 내 유전자의 mRNA 양이 많다고 해서 반드시 단백질이 많이 만들어지는 것은 아닙니다. `Lin28a`가 억제(siLin28a)되었을 때, 전사체 수준(RNA-seq)의 변화와 실제 리보솜이 번역 중인 수준(RPF)의 변화를 비교하여 **유전자별 번역 효율(TE)의 증감**을 정량적으로 규명합니다.

#### 2) 상세 Method (수행 방법)

1. **발현량 행렬 도출**: 1주차에서 사용한 `featureCounts`를 확장 적용하여 `RNA-control(siLuc)`, `RNA-siLin28a`, `RPF-control(siLuc)`, `RPF-siLin28a` 등 모든 조건의 유전자별 Read Count 테이블을 생성합니다.
2. **데이터 정규화 (Normalization)**: 총 Read 수의 차이를 보정하기 위해 파이썬(`pandas`)을 이용해 각 샘플의 카운트를 **TPM(Transcripts Per Million)** 또는 **CPM(Counts Per Million)** 값으로 변환합니다.
3. **번역 효율(TE) 계산**: 각 유전자마다 아래 식을 이용해 번역 효율을 구합니다.

$$\text{TE} = \frac{\text{RPF 정규화 값}}{\text{RNA 정규화 값}}$$


4. **상관관계 및 시각화**:
* X축을 $\log_2(\text{RNA의 변화량})$, Y축을 $\log_2(\text{RPF의 변화량})$으로 하는 Scatter plot을 그립니다.
* 컨트롤 대비 siLin28a에서 TE가 유의미하게 증가하거나 감소한 유전자 군(Outlier)을 동정합니다.


5. **인사이트 확장**: 1주차에 결합했던 **Cellular Localization 데이터**와 연계하여, "Lin28a 넉다운 시 주로 세포질(Cytoplasm)에 있는 유전자들의 번역 효율이 저하되는가?"와 같은 공간적 특성을 통계적으로 검증합니다.

---

### 🚀 추천 2: 리보솜 주기성(Ribosome Periodicity) 분석 및 Triplet Phasing 확인

* **연계 실습**: 2주차(bedtools를 이용한 개시 코돈 주변 5' 말단 카운트 가공)

#### 2) 실습 목적

리보솜은 mRNA 위를 임의로 움직이지 않고, 번역틀(Reading Frame)을 유지하며 **3개 염기(1코돈)씩 불연속적으로 이동**합니다. 이 생물학적 현상(Triplet Phasing)이 실제 우리가 정제한 RPF(Ribosome Protected Fragment) 데이터의 5' 말단 분포에서도 3bp 주기의 피크로 나타나는지 데이터로 직접 증명합니다.

#### 2) 상세 Method (수행 방법)

1. **상대 좌표 변환**: 2주차 실습에서 `bedtools intersect`를 통해 얻은 `fivepcounts-filtered-RPF-siLuc.txt` 데이터를 파이썬으로 로드합니다. 유전체 절대 좌표를 **개시 코돈(Start Codon)의 첫 번째 염기 자리를 `0`으로 하는 상대 좌표**로 변환합니다.
2. **리드 카운트 합산 (Aggregation)**: 모든 유전자에 대해 개시 코돈 주변(예: -30bp ~ +90bp 영역)의 상대적 위치별 리드 카운트를 전부 더해 하나의 테이블로 메타화합니다.
3. **메타 유전자 플롯(Meta-gene plot) 시각화**: X축을 개시 코돈으로부터의 거리(Base), Y축을 합산된 리드 수로 하여 바 차트(Bar chart) 또는 라인 플롯을 그립니다.
4. **주기성 통계 검증**:
* 개시 코돈 이후(+0, +3, +6, +9...) 3의 배수 자리에서 리드가 주기적으로 폭발하는지 확인합니다.
* 파이썬의 `scipy.signal` 패키지를 활용해 Autocorrelation(자기상관분석)을 수행하거나 FFT(고속 푸리에 변환)를 적용하여 데이터 내에 '주기 3'의 신호가 통계적으로 확실히 존재하는지 분석합니다.



---

### 🚀 추천 3: CLIP-seq 데이터를 활용한 Lin28a 결합 부위의 Shannon Entropy 및 모티프(Motif) 발굴

* **연계 실습**: 3주차(염기별 카운팅 및 Shannon Entropy 계산, UCSC 브라우저 시각화)

#### 1) 실습 목적

3주차 실습에서는 특정 타겟 유전자 구간의 Shannon Entropy를 구했습니다. CLIP-seq 데이터는 단백질이 RNA에 결합한 부위를 캡처하는데, 실험 과정(UV Cross-linking 및 역전사)의 특성상 **단백질이 결합했던 염기 자리에 돌연변이(Mutation)나 결실(Deletion)이 집중**되어 Shannon Entropy가 비정상적으로 높아집니다. 이 성질을 이용해 Lin28a가 실제 결합하는 RNA 상의 핵심 서열(Motif)을 찾아냅니다.

#### 2) 상세 Method (수행 방법)

1. **전장 혹은 특정 염색체 엔트로피 스캔**: 3주차에 작성한 파이썬 코드를 확장하여, `CLIP-35L33G.bam` 파일로부터 리드가 많이 매핑된 주요 유전자 영역 전체의 단일 염기별 Shannon Entropy를 자동 계산합니다.
2. **엔트로피 피크(Cross-linking Site) 검출**: Shannon Entropy 값이 특정 임계값(Threshold) 이상으로 높게 튀는 고위험(?) 변이 좌표들을 파이썬 코드로 필터링하여 추출합니다. 이 좌표들이 곧 Lin28a가 강하게 결합했던 자리입니다.
3. **주변 서열(Context Window) 추출**: 추출된 변이 좌표를 중심으로 좌우 10~20bp에 해당하는 유전체 서열(Reference Sequence, `mm39.fa`)을 `bedtools getfasta` 또는 파이썬 바이오파이썬(Biopython) 툴을 이용해 FASTA 파일로 일괄 추출합니다.
4. **서열 모티프 분석**:
* 추출한 서열들을 웹 기반 모티프 발굴 툴인 MEME Suite(DREME/MEME)에 업로드합니다.
* 또는 파이썬으로 직접 4-mer나 5-mer 서열 빈도를 계산하여, 엔트로피가 높은 곳 주변에 유독 자주 등장하는 서열 패턴(예: `GGAGA` 등 기존 연구에 알려진 Lin28a 결합 모티프)이 존재하는지 비교·검증하고 이를 WebLogo 형태로 시각화합니다.